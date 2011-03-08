package ibis.amuse;

import ibis.util.ThreadPool;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.IOException;
import java.net.Socket;
import java.nio.ByteBuffer;
import java.nio.channels.SocketChannel;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class Worker implements Runnable {

	private static final Logger logger = LoggerFactory.getLogger(Worker.class);

	private final SocketChannel channel;
	
	Worker(SocketChannel socket) throws IOException {
		this.channel = socket;

		logger.info("New connection from " + socket.socket().getRemoteSocketAddress());

		ThreadPool.createNew(this, "Amuse worker");
	}
	
	private void initialize() throws IOException {
		
		Call init = new Call(channel);
//		
//		String magic = readString();
//		logger.debug("magic string = " + magic);
//
//		String name = readString();
//		logger.debug("worker name = " + name);
//
//		String hostname = readString();
//		logger.debug("hostname = " + hostname);
//
//		writeString("OK");
//		out.flush();
//
//		try {
//			Thread.sleep(1000);
//		} catch (InterruptedException e) {
//			// IGNORE
//		}
//		writeString("STARTED");
//		out.flush();
	}

	public void run() {
		try {
			initialize();
		} catch (IOException e) {
			logger.error("Exception in worker", e);
		}

	}

}
