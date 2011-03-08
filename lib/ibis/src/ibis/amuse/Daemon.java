package ibis.amuse;

import java.io.IOException;
import java.net.InetSocketAddress;
import java.net.ServerSocket;
import java.net.Socket;
import java.net.SocketAddress;
import java.nio.channels.ServerSocketChannel;
import java.nio.channels.SocketChannel;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Ibis Daemon. Starts workers for Amuse using software from the Ibis project.
 * 
 * @author Niels Drost
 * 
 */
public class Daemon {

	private static final Logger logger = LoggerFactory.getLogger(Daemon.class);

	public static void main(String[] arguments) throws IOException {
		ServerSocketChannel serverSocket = ServerSocketChannel.open();
		serverSocket.socket().bind(new InetSocketAddress(61575));

		while (true) {
			try {
				logger.info("Waiting for connection");
				SocketChannel socket = serverSocket.accept();

				new Worker(socket);
			} catch (Exception e) {
				logger.error("Error while handling new worker connection", e);

				// wait a bit before handling the next connection
				try {
					Thread.sleep(1000);
				} catch (InterruptedException e1) {
					// IGNORE
				}
			}
		}

	}

}
