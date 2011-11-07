package ibis.amuse;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;

public class OutputPrefixForwarder extends Thread {

    private final BufferedReader input;

    private final PrintStream output;

    private final String outputPrefix;

    /**
     * @param input
     *            Input stream
     * @param output
     *            Stream to forward output to
     * @param outputPrefix
     *            Prefix to add to all lines of output
     * 
     * @throws IOException
     *             if the reading stream cannot be created.
     */
    public OutputPrefixForwarder(InputStream input, PrintStream output,
            String outputPrefix) throws IOException {
        this.input = new BufferedReader(new InputStreamReader(input));
        this.output = output;

        this.outputPrefix = outputPrefix;

        setDaemon(false);
        setName("prefix forwarder");
        start();
    }

    /**
     * Forwards input stream to given output stream.
     */
    public void run() {
        while (true) {
            try {
                String line = input.readLine();

                if (line == null) {
                    // we're done
                	output.flush();
                    return;
                }

                output.println(outputPrefix + line);
            } catch (IOException e) {
                // we're done
            	output.flush();
                return;
            }
        }

    }

}
