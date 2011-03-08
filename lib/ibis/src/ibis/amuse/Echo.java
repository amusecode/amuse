package ibis.amuse;
import java.io.IOException;
import java.net.ServerSocket;
import java.net.Socket;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class Echo {
	
	private static final Logger logger = LoggerFactory.getLogger(Echo.class);

	public static void handleConnection(Socket socket) {

		byte[] buffer = new byte[100];

		try {
			while (true) {
				logger.info("Waiting for input");
				int read = socket.getInputStream().read(buffer);

				if (read == -1) {
					logger.info("Done!");
					return;
				}

				for (int i = 0; i < read; i++) {
					System.err.format("%x\n", buffer[i]);
				}

				socket.getOutputStream().write(buffer, 0, read);
			}

		} catch (IOException e) {
			logger.error("Error on handling connection", e);
		} finally {
			try {
				socket.close();
			} catch (IOException e) {
				//IGNORE
			}
		}

	}

	public static void main(String[] arguments) throws IOException {
		ServerSocket serverSocket = new ServerSocket(61575);

		while (true) {

			logger.info("Waiting for connection");
			Socket socket = serverSocket.accept();

			logger.info("New connection from " + socket.getRemoteSocketAddress());
			
			handleConnection(socket);
		}

	}

}
