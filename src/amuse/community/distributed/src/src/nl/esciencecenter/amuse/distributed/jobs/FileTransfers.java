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
package nl.esciencecenter.amuse.distributed.jobs;

import ibis.ipl.ReadMessage;
import ibis.ipl.WriteMessage;

import java.io.File;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Utility class for transferring files, mostly to/from ipl messages
 * 
 * @author Niels Drost
 * 
 */
public class FileTransfers {
    
    public static final int BUFFER_SIZE = 1000;

    private static final Logger logger = LoggerFactory.getLogger(FileTransfers.class);

    static void writeFile(File file, WriteMessage writeMessage, ByteBuffer buffer) throws IOException {
        writeMessage.writeString(file.getPath());
        writeMessage.writeLong(file.length());
        try (FileChannel channel = FileChannel.open(file.toPath(), StandardOpenOption.READ)) {
            while (true) {
                //write file in 1000 byte chunks in case it is big
                buffer.clear();

                int read = channel.read(buffer);

                if (read == -1) {
                    return;
                }

                buffer.flip();
                writeMessage.writeByteBuffer(buffer);
            }
        }

    }

    static void addFiles(File file, ArrayList<File> result) {
        if (file.isFile()) {
            result.add(file);
            logger.debug("Added file: " + file);
        } else if (file.isDirectory()) {
            for (File child : file.listFiles()) {
                addFiles(child, result);
            }
        }
    }

    static void writeDirectory(String filename, WriteMessage writeMessage) throws IOException {
        ByteBuffer buffer = ByteBuffer.allocate(BUFFER_SIZE);
        
        File directory = new File(filename);

        if (!directory.isDirectory()) {
            throw new IOException("Directory \"" + filename + "\" not found");
        }

        ArrayList<File> files = new ArrayList<File>();

        //recursively add all files
        addFiles(directory, files);

        writeMessage.writeInt(files.size());
        for (File file : files) {
            writeFile(file, writeMessage, buffer);
        }
    }

    static void readFile(File directory, ReadMessage readMessage, ByteBuffer buffer) throws IOException {
        String filename = readMessage.readString();
        long size = readMessage.readLong();

        File file = new File(directory, filename);

        logger.debug("Reading file " + file);

        try (FileChannel channel = FileChannel.open(file.toPath(), StandardOpenOption.WRITE, StandardOpenOption.CREATE)) {

            long bytesLeft = size;

            while (bytesLeft > 0) {
                buffer.clear();
                buffer.limit((int) Math.min(buffer.capacity(), bytesLeft));

                readMessage.readByteBuffer(buffer);

                buffer.flip();

                while (buffer.hasRemaining()) {
                    int written = channel.write(buffer);
                    bytesLeft -= written;
                }
            }
        }
    }

    static void readDirectory(String filename, ReadMessage readMessage) throws IOException {
        ByteBuffer buffer = ByteBuffer.allocate(BUFFER_SIZE);
        
        File directory = new File(filename);

        if (!directory.isDirectory()) {
            directory.mkdir();
        }

        if (!directory.isDirectory()) {
            throw new IOException("Directory \"" + filename + "\" could not be created");
        }

        int count = readMessage.readInt();

        for (int i = 0; i < count; i++) {
            readFile(directory, readMessage, buffer);
        }

    }

}
