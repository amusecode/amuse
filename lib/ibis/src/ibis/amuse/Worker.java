package ibis.amuse;

import ibis.util.ThreadPool;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.channels.SocketChannel;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class Worker implements Runnable {

    private static final Logger logger = LoggerFactory.getLogger(Worker.class);

    private final SocketChannel channel;

    private String name;

    private String hostname;

    Worker(SocketChannel socket) throws IOException {
        this.channel = socket;

        logger.info("New connection from "
                + socket.socket().getRemoteSocketAddress());

        ThreadPool.createNew(this, "Amuse worker");
    }

    public synchronized String getName() {
        return name;
    }

    public synchronized String getHostname() {
        return hostname;
    }

    private void initialize() throws IOException {

        ByteBuffer magic = ByteBuffer.allocate(MessageFormat.MAGIC_STRING
                .length());
        Call.readAll(channel, magic);
        String receivedString = new String(magic.array(), "UTF-8");
        if (!receivedString.equals(MessageFormat.MAGIC_STRING)) {
            throw new IOException("magic string " + MessageFormat.MAGIC_STRING
                    + " not received. Instead got: " + receivedString);
        }

        Call init = new Call(channel);

        synchronized (this) {
            this.name = init.getString(0);
            this.hostname = init.getString(1);
        }

        logger.info("Initializing " + this);

        Result result = new Result(init.getId(), init.getFunctionID(),
                init.getCount());

        result.writeTo(channel);
    }

    public void run() {
        try {
            initialize();
        } catch (IOException e) {
            logger.error("Exception in worker", e);
        }

        while (true) {
            Call call;
            try {
                call = new Call(channel);

                Result result = null;
                if (call.getFunctionID() == MessageFormat.FUNCTION_ID_STOP) {
                    result = end(call);
                    
                    //write result, close channel, end thread
                    result.writeTo(channel);
                    channel.close();
                    logger.info(this + " done!");
                    return;
                } else if (call.getFunctionID() == MessageFormat.FUNCTION_ID_REDIRECT_OUTPUT) {
                    result = setLogFiles(call);
                } else {
                    result = new Result(call.getId(), call.getFunctionID(),
                            new Exception("cannot actually do any calls yet!"));
                }
                result.writeTo(channel);
            } catch (IOException e) {
                logger.error("Error on handling call, stopping worker", e);
                return;
            }

        }

    }

    private Result end(Call call) {
        Result result = new Result(call.getId(), call.getFunctionID(),
                call.getCount());
        return result;
    }

    private Result setLogFiles(Call call) throws IOException {
        try {
            logger.info("Setting stdout to " + call.getString(0)
                    + ", and stderr to " + call.getString(1));
            Result result = new Result(call.getId(), call.getFunctionID(),
                    call.getCount(), new int[1], new long[0], new float[0],
                    new double[0], new boolean[0], new String[0]);
            return result;
        } catch (Throwable t) {
            // return exception as result
            Result result = new Result(call.getId(), call.getFunctionID(), t);
            return result;
        }
    }

    public String toString() {
        return "Worker \"" + getName() + "@" + getHostname() + "\"";
    }
}
