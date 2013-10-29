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
package nl.esciencecenter.amuse.distributed.dummycode;

import java.io.IOException;
import java.net.InetSocketAddress;
import java.net.SocketAddress;
import java.nio.ByteBuffer;
import java.nio.channels.SocketChannel;
import java.nio.charset.StandardCharsets;

import nl.esciencecenter.amuse.distributed.AmuseMessage;

/**
 * Class that mimics an amuse script starting a worker and making calls. Only uses the "Dummy" code included in distributed amuse.
 * Useful mostly for testing.
 * 
 * @author Niels Drost
 * 
 */
public class FakeAmuse {

    public static final String MAGIC_STRING = "TYPE_WORKER";
    
    private final SocketChannel channel;
    
    // Fakes the following python code:
    //        self.socket.connect((self.daemon_host, self.daemon_port))
    //        
    //        self.socket.setblocking(1)
    //        
    //        self.socket.setsockopt(socket.SOL_TCP, socket.TCP_NODELAY, 1)
    //        
    //        self.socket.sendall('TYPE_WORKER'.encode('utf-8'))
    //        
    //        arguments = {'string': [self.name_of_the_worker, self.worker_dir, self.hostname, self.redirect_stdout_file, 
    //        self.redirect_stderr_file, self.node_label], 
    //        'int32': [self.number_of_workers, self.number_of_nodes, self.number_of_threads], 
    //        'bool': [self.copy_worker_code]}
    //        
    //        message = SocketMessage(call_id=1, function_id=10101010, call_count=1, dtype_to_arguments=arguments);
    //
    //        message.send(self.socket)
    //        
    //        logging.getLogger("channel").info("waiting for worker %s to be initialized", self.name_of_the_worker)
    //
    //        result = SocketMessage()
    //        result.receive(self.socket)
    //        
    //        if result.error:
    //            logging.getLogger("channel").error("Could not start worker: %s", result.strings[0])
    //            self.stop()
    //            raise exceptions.CodeException("Could not start worker for " + self.name_of_the_worker + ": " + result.strings[0])
    //        
    //        logging.getLogger("channel").info("worker %s initialized", self.name_of_the_worker)
    public FakeAmuse(int port) throws IOException {

        SocketAddress address = new InetSocketAddress("localhost", port);
        channel = SocketChannel.open(address);

        byte[] magicBytes = MAGIC_STRING.getBytes(StandardCharsets.UTF_8);
        
        ByteBuffer magicBuffer = ByteBuffer.wrap(magicBytes);
        
        channel.write(magicBuffer);
        
        AmuseMessage request = new AmuseMessage(1,  10101010, 1);
        
        //self.name_of_the_worker, self.worker_dir, self.hostname, self.redirect_stdout_file, self.redirect_stderr_file, self.node_label
        
        String[] strings = { "fake", "fake", "localhost", "/dev/null", "/dev/null", "default"};
        
        
    }
    
    
    
}
