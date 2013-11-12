/*
 * Copyright 2013 Netherlands eScience Center
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package nl.esciencecenter.amuse.distributed.web;

import java.io.IOException;
import java.io.PrintWriter;
import java.net.InetAddress;
import java.net.UnknownHostException;
import java.util.Map;

import javax.servlet.ServletException;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import nl.esciencecenter.amuse.distributed.DistributedAmuse;
import nl.esciencecenter.amuse.distributed.DistributedAmuseException;
import nl.esciencecenter.amuse.distributed.jobs.Job;
import nl.esciencecenter.amuse.distributed.jobs.WorkerJob;
import nl.esciencecenter.amuse.distributed.reservations.Reservation;
import nl.esciencecenter.amuse.distributed.resources.Resource;

import org.eclipse.jetty.server.Request;
import org.eclipse.jetty.server.Server;
import org.eclipse.jetty.server.ServerConnector;
import org.eclipse.jetty.server.handler.AbstractHandler;
import org.eclipse.jetty.util.thread.QueuedThreadPool;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * @author Niels Drost
 * 
 *         Small web interface to monitor distributed amuse.
 * 
 */
public class WebInterface extends AbstractHandler {

    private static final Logger logger = LoggerFactory.getLogger(WebInterface.class);

    private final int port;

    private final Server server;

    private final DistributedAmuse distributedAmuse;

    private static void writeMapAsTable(Map<String, String> map, PrintWriter writer) {
        writer.println("<table cellspacing=\"10\">");

        for (Map.Entry<String, String> entry : map.entrySet()) {
            writer.println("<tr><td><pre>" + entry.getKey() + "</pre></td><td><pre>" + entry.getValue() + "</pre></td></tr>");
        }

        writer.println("</table>");
    }

    public WebInterface(DistributedAmuse distributedAmuse, int port) throws Exception {
        this.distributedAmuse = distributedAmuse;

        QueuedThreadPool threadPool = new QueuedThreadPool(100);
        threadPool.setDaemon(true);
        
        server = new Server(threadPool);
        
        
        ServerConnector http = new ServerConnector(server);
        //http.setHost("localhost");
        http.setPort(port);
        //http.setIdleTimeout(30000);
        
        server.addConnector(http);
        server.setHandler(this);

        server.start();

        //get actual port from jetty
        this.port = ((ServerConnector) server.getConnectors()[0]).getLocalPort();
        logger.info("Web interface running on " + getURL());
    }

    public int getPort() {
        return port;
    }
    
    public String getURL() {
        try {
            InetAddress localAddress = InetAddress.getLocalHost();
            return "http://" + localAddress.getHostAddress() + ":" + this.port;
        } catch (UnknownHostException e) {
            return "http://localhost:" + this.port;
        }
    }

    @Override
    public void handle(String target, Request baseRequest, HttpServletRequest request, HttpServletResponse response)
            throws IOException, ServletException {
        logger.debug("handling request for page: " + target);
        response.setContentType("text/html;charset=utf-8");
        writeHeader(response.getWriter());
        try {
            if (target.equals("/")) {
                writeOverview(response.getWriter());
            } else if (target.equals("/resources")) {
                writeResourceTable(response.getWriter());
            } else if (target.equals("/reservations")) {
                writeReservationTable(response.getWriter());
            } else if (target.equals("/jobs")) {
                writeJobTable(response.getWriter());
            } else if (target.startsWith("/resources/")) {
                String resourceID = target.substring("/resources/".length());
                writeResourceDetailsResponse(response.getWriter(), resourceID);
            } else if (target.startsWith("/reservations/")) {
                String reservationID = target.substring("/reservations/".length());
                writeReservationDetailsResponse(response.getWriter(), reservationID);
            } else if (target.startsWith("/jobs/")) {
                String jobID = target.substring("/jobs/".length());
                writeJobDetailsResponse(response.getWriter(), jobID);
            } else {
                throw new DistributedAmuseException("unknown web interface path: " + target);
            }
            response.setStatus(HttpServletResponse.SC_OK);
        } catch (DistributedAmuseException e) {
            response.getWriter().println("<pre>");
            e.printStackTrace(response.getWriter());
            response.getWriter().println("</pre>");
            response.setStatus(HttpServletResponse.SC_INTERNAL_SERVER_ERROR);
        }
        baseRequest.setHandled(true);
        writeFooter(response.getWriter());
    }

    private static void writeHeader(PrintWriter writer) {
        writer.println("<html>");
        writer.println("<head>");
        writer.println("<title>Distributed Amuse Web Interface</title>");
        writer.println("</head>");
        writer.println("<body>");
    }

    private static void writeFooter(PrintWriter writer) {
        writer.println("</body>");
        writer.println("</html>");
    }

    private void writeOverview(PrintWriter writer) throws DistributedAmuseException {
        writer.println("<h1>Distributed Amuse Webinterface</h1>");

        writer.println("<h2>Resources</h2>");
        writeResourceTable(writer);

        writer.println("<h2>Reservations</h2>");
        writeReservationTable(writer);

        writer.println("<h2>Jobs</h2>");
        writeJobTable(writer);

    }

    private void writeResourceTable(PrintWriter writer) {
        writer.println("<table border=\"1px\">");

        writer.println("<tr><th>ID</th><th>Name</th><th>Hostname</th><th>Type</th></tr>");
        for (Resource resource : distributedAmuse.resourceManager().getResources()) {
            writer.printf("<tr><td><a href=/resources/%s>%s</a></td><td>%s</td><td>%s</td><td>%s</td></tr>\n", resource.getId(),
                    resource.getId(), resource.getName(), resource.getLocation(), resource.getSchedulerType());
        }
        writer.println("</table>");
    }

    private void writeReservationTable(PrintWriter writer) throws DistributedAmuseException {
        writer.println("<table border=\"1px\">");

        writer.println("<tr><th>ID</th><th>Node Label</th><th>Resource Name</th><th>Node Count</th><th>Status</th></tr>");
        for (Reservation reservation : distributedAmuse.reservationManager().getReservations()) {
            String state = distributedAmuse.reservationManager().getState(reservation);
            writer.printf("<tr><td><a href=/reservations/%s>%s</a></td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>\n",
                    reservation.getID(), reservation.getID(), reservation.getNodeLabel(), reservation.getResourceName(),
                    reservation.getNodeCount(), state);
        }
        writer.println("</table>");
    }

    private void writeJobTable(PrintWriter writer) {
        writer.println("<table border=\"1px\">");

        writer.println("<tr><th>ID</th><th>Job State</th><th>Job Type</th><th>Job Label</th></tr>");
        for (WorkerJob job : distributedAmuse.jobManager().getWorkerJobs()) {
            writer.printf("<tr><td><a href=/jobs/%d>%d</a></td><td>%s</td><td>%s</td><td>%s</td></tr>\n", job.getJobID(),
                    job.getJobID(), job.getJobState(), "worker", job.getNodeLabel());
        }
        writer.println("</table>");
    }

    private void writeResourceDetailsResponse(PrintWriter writer, String resourceID) throws DistributedAmuseException {
        int id = Integer.parseInt(resourceID);

        Resource resource = distributedAmuse.resourceManager().getResource(id);

        writer.println("<h1>Resource " + resourceID + "</h1>");
        writeMapAsTable(resource.getStatusMap(), writer);
    }

    private void writeReservationDetailsResponse(PrintWriter writer, String reservationID) throws DistributedAmuseException {
        int id = Integer.parseInt(reservationID);

        Reservation reservation = distributedAmuse.reservationManager().getReservation(id);

        writer.println("<h1>Reservation " + reservationID + "</h1>");

        Map<String, String> statusMap = reservation.getStatusMap();

        String resourceID = statusMap.get("Resource");

        //make resource value into a link
        statusMap.put("Resource", "<a href=../resources/" + resourceID + ">" + resourceID + "</a>");

        writeMapAsTable(statusMap, writer);

        writer.println("<h3> Standard output</h3>");
        writer.println("<p>");
        try {
            for (String line : reservation.getStdout()) {
                writer.println(line + "<br>");
            }
        } catch (DistributedAmuseException e) {
            writer.println("Failed to get standard output: " + e.getMessage());
            writer.println("<pre>");
            e.printStackTrace(writer);
            writer.println("</pre>");

        }
        writer.println("</p>");

        writer.println("<h3> Standard error</h3>");
        writer.println("<p>");
        try {
            for (String line : reservation.getStderr()) {
                writer.println(line + "<br>");
            }
        } catch (DistributedAmuseException e) {
            writer.println("Failed to get standard error: " + e.getMessage());
            writer.println("<pre>");
            e.printStackTrace(writer);
            writer.println("</pre>");
        }

        writer.println("</p>");

    }

    private void writeJobDetailsResponse(PrintWriter writer, String jobID) throws DistributedAmuseException {
        int id = Integer.parseInt(jobID);

        Job job = distributedAmuse.jobManager().getJob(id);
        writer.println("<h1>Job " + jobID + "</h1>");
        writeMapAsTable(job.getStatusMap(), writer);
    }

    /**
     * Stop web interface
     */
    public void end() {
        try {
            server.stop();
        } catch (Exception e) {
            logger.error("Error while stopping web interface", e);
        }
    }
}
