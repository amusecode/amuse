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
package nl.esciencecenter.amuse.distributed.util;

import ibis.ipl.ReadMessage;
import ibis.ipl.WriteMessage;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;

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

    private static void writeFile(Path file, Path root, WriteMessage writeMessage, ByteBuffer buffer) throws IOException {
        String relativePath = root.relativize(file).toString();
        long size = Files.size(file);

        logger.debug("Writing file {} with size {} and root {} to relative path {}", file, size, root, relativePath);

        //not a directory
        writeMessage.writeBoolean(false);
        writeMessage.writeString(relativePath);
        writeMessage.writeLong(size);
        try (FileChannel channel = FileChannel.open(file, StandardOpenOption.READ)) {
            while (true) {
                //write file in 1000 byte chunks in case it is big
                buffer.clear();

                int read = channel.read(buffer);

                if (read == -1) {
                    return;
                }

                buffer.flip();
                writeMessage.writeArray(buffer.array(), 0, read);
                writeMessage.flush();
            }
        }

    }

    private static void writeDirectory(Path file, Path root, WriteMessage writeMessage) throws IOException {
        String relativePath = root.relativize(file).toString();

        logger.debug("Writing directory entry with path {} and root {} to relative path {}", file, root, relativePath);

        //true: this is a directory
        writeMessage.writeBoolean(true);
        writeMessage.writeString(relativePath);
        writeMessage.writeLong(0);

    }

    private static void writeItem(Path file, Path root, WriteMessage writeMessage, ByteBuffer buffer) throws IOException {
        logger.debug("Writing item {} with root {}", file, root);
        if (Files.isDirectory(file)) {
            writeDirectory(file, root, writeMessage);

            for (Path child : Files.newDirectoryStream(file)) {
                writeItem(child, root, writeMessage, buffer);
            }
        } else if (Files.isRegularFile(file)) {
            writeFile(file, root, writeMessage, buffer);
        }
    }

    public static void writeDirectory(String directoryName, Path root, WriteMessage writeMessage) throws IOException {
        logger.debug("Writing directory with name {} and root {}", directoryName, root);

        ByteBuffer buffer = ByteBuffer.allocate(BUFFER_SIZE);

        Path directory = root.resolve(directoryName);

        if (!Files.isDirectory(directory)) {
            throw new IOException("Directory \"" + directory + "\" not found");
        }

        writeItem(directory, root, writeMessage, buffer);

        //signal this was the last file
        writeMessage.writeBoolean(false);
        writeMessage.writeString("end");
        writeMessage.writeLong(-1);
    }

    public static void writeFilesInDirectory(Path directory, WriteMessage writeMessage, String filePattern) throws IOException {
        logger.debug("Writing directory {}", directory);

        ByteBuffer buffer = ByteBuffer.allocate(BUFFER_SIZE);

        if (!Files.isDirectory(directory)) {
            throw new IOException("Directory \"" + directory + "\" not found");
        }

        //write all files in this directory
        for (Path child : Files.newDirectoryStream(directory)) {
            if (Files.isRegularFile(child) && child.getFileName().toString().matches(filePattern)) {
                writeFile(child, directory, writeMessage, buffer);
            }
        }

        //signal this was the last file
        writeMessage.writeBoolean(false);
        writeMessage.writeString("end");
        writeMessage.writeLong(-1);
    }

    private static void readFile(Path path, long size, ReadMessage readMessage, ByteBuffer buffer) throws IOException {
        logger.debug("Reading file with path {} and size {}", path, size);

        try (FileChannel channel = FileChannel.open(path, StandardOpenOption.WRITE, StandardOpenOption.CREATE)) {

            long bytesLeft = size;

            while (bytesLeft > 0) {
                int bytesToRead = (int) Math.min(buffer.capacity(), bytesLeft);

                readMessage.readArray(buffer.array(), 0, bytesToRead);

                buffer.clear();
                buffer.limit(bytesToRead);
                
                while (buffer.hasRemaining()) {
                    int written = channel.write(buffer);
                    bytesLeft -= written;
                }
            }
        }
    }

    public static void readDirectory(Path root, ReadMessage readMessage) throws IOException {
        logger.debug("Reading directory contents with root {}", root);

        ByteBuffer buffer = ByteBuffer.allocate(BUFFER_SIZE);

        if (!Files.isDirectory(root)) {
            throw new IOException("Directory \"" + root + "\" does not exist");
        }

        while (true) {
            boolean isDirectory = readMessage.readBoolean();
            String filename = readMessage.readString();
            long size = readMessage.readLong();

            logger.debug("Reading item {} with filename {} and size {}", isDirectory ? "directory" : "file", filename, size);

            if (filename.equals("end") && size == -1 && isDirectory == false) {
                //last file marker
                logger.debug("Reading done");
                return;
            }

            Path path = root.resolve(filename);

            if (isDirectory) {
                logger.debug("Creating directory {}", path);
                if (!Files.isDirectory(path)) {
                    Files.createDirectory(path);
                }
            } else {
                readFile(path, size, readMessage, buffer);
            }

        }

    }

}
