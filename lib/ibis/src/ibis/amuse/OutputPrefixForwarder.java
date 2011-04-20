package ibis.amuse;

import ibis.util.ThreadPool;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;

public class OutputPrefixForwarder implements Runnable {

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

        ThreadPool.createNew(this, "prefix forwarder");
    }

    /**
     * Forwards standard out of server to given output stream. Filters out line
     * containing server address.
     */
    public void run() {
        while (true) {
            try {
                String line = input.readLine();

                if (line == null) {
                    // we're done
                    return;
                }

                output.println(outputPrefix + line);
            } catch (IOException e) {
                // we're done
                return;
            }
        }

    }

}
