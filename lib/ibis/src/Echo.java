import java.io.IOException;
import java.net.ServerSocket;
import java.net.Socket;

public class Echo {

    public static void main(String[] arguments) throws IOException {
        ServerSocket serverSocket = new ServerSocket(61575);

        System.err.println("Waiting for connection");
        Socket socket = serverSocket.accept();

        byte[] buffer = new byte[100];

        while (true) {
            System.err.println("Waiting for input");
            int read = socket.getInputStream().read(buffer);

            if (read == -1) {
                System.err.println("Done!");
                return;
            }

            for (int i = 0; i < read; i++) {
                System.err.format("%x\n", buffer[i]);
            }

            socket.getOutputStream().write(buffer, 0, read);
        }

    }

}
