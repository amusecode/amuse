package nl.esciencecenter.amuse.distributed.dummycode;


public class Code implements CodeInterface {

  public Code(String codeDir, String amuseRootDir) {
  }

  public void end() {
        //IGNORE
  }
  
  // parameter "int_in" is an input parameter
  // parameter "int_out" is an output parameter
  // parameter "len" is a length parameter
  public int echo_array_with_result(int[] int_in, int[] int_out, int len) {
	for(int i = 0; i < len; i++) {
	    int_out[i] = int_in[i];
	}
	return -1;
  }
  
  
  // parameter "string_in" is an input parameter
  public int print_error_string(String[] string_in) {
	for(int i = 0; i < string_in.length; i++) {
    	    System.err.println(string_in[i]);
	}
	return 0;
  }
  
  
  // parameter "string_inout1" is an inout parameter
  // parameter "string_inout2" is an inout parameter
  public int echo_strings(String[] string_inout1, String[] string_inout2) {
	for(int i = 0; i < string_inout1.length; i++) {
	    String tmp = string_inout1[i];
	    string_inout1[i] = string_inout2[i];
	    string_inout2[i] = tmp;
	}
	return 0;
  }
  
  
  // parameter "double_in" is an input parameter
  // parameter "double_out" is an output parameter
  public int echo_double(double[] double_in, double[] double_out) {
	for(int i = 0; i < double_in.length; i++) {
	    double_out[i] = double_in[i];
	}
	return 0;
  }
  
  
  // parameter "float_in" is an input parameter
  // parameter "float_out" is an output parameter
  public int echo_float(float[] float_in, float[] float_out) {
	for(int i = 0; i < float_in.length; i++) {
	    float_out[i] = float_in[i];
	}
	return 0;
  }
  
  
  // parameter "int_in" is an input parameter
  // parameter "int_out" is an output parameter
  public int echo_int(int int_in, int[] int_out) {
	int_out[0] = int_in;
	if (int_in < 0) {
	    return -1;
	} else {
   	    return 0;
	}
  }

  // parameter "int_in" is an input parameter
  // parameter "int_out" is an output parameter
  public int echo_int_array(int[] int_in, int[] int_out) {
	for(int i = 0; i < int_in.length; i++) {
	    int_out[i] = int_in[i];
	}
	return 0;
  }
   
  
  // parameter "in_out" is an inout parameter
  // parameter "len" is a length parameter
  public int echo_inout_array_with_result(int[] in_out, int len) {
      for(int i = 0; i < len; i++) {
          in_out[i] = in_out[i] + 10;
      }
      return 11;
  }
  
  
  // parameter "int_in" is an input parameter
  // parameter "int_out" is an output parameter
  // parameter "len" is a length parameter
  public void echo_array(int[] int_in, int[] int_out, int len) {
	for(int i = 0; i < len; i++) {
	    int_out[i] = int_in[i];
	}
  }
  
  
  // parameter "in" is an input parameter
  // parameter "out" is an output parameter
  public int echo_long_long_int(long[] in, long[] out) {
	for(int i = 0; i < in.length; i++) {
 	    out[i] = in[i];
	}
	if (in[0] < 0) {
	    return -1;
	} else {
   	    return 0;
	}
 }
  
  
  // parameter "string_in" is an input parameter
  // parameter "string_out" is an output parameter
  public int echo_string(String[] string_in, String[] string_out) {
	for(int i = 0; i < string_in.length; i++) {
	    string_out[i] = string_in[i];
	}
	return 0;
  }
  
  
  // parameter "input" is an input parameter
  // parameter "output" is an output parameter
  public int echo_logical(boolean[] input, boolean[] output) {
	for(int i = 0; i < input.length; i++) {
	    output[i] = input[i];
	}
	return 0;
  }
  
  
  // parameter "string_in" is an input parameter
  public int print_string(String[] string_in) {
	for(int i = 0; i < string_in.length; i++) {
    	    System.out.println(string_in[i]);
	}
	return 0;
  }
  
  
  // parameter "int_in1" is an input parameter
  // parameter "int_in2" is an input parameter
  // parameter "int_out1" is an output parameter
  // parameter "int_out2" is an output parameter
  // parameter "len" is a length parameter
  public int echo_2_int(int[] int_in1, int[] int_in2, int[] int_out1, int[] int_out2, int len) {
	for(int i = 0; i < len; i++) {
            int_out1[i] = int_in1[i];
            int_out2[i] = int_in2[i];
    	}

	return len;
  }
  
}
