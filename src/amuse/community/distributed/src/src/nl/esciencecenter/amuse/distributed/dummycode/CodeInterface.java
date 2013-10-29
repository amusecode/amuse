package nl.esciencecenter.amuse.distributed.dummycode;

public interface CodeInterface {
  public void end();
  
  
  // parameter "int_in" is an input parameter
  // parameter "int_out" is an output parameter
  int echo_int(int int_in, int[] int_out);
  
  
  // parameter "int_in" is an input parameter
  // parameter "int_out" is an output parameter
  // parameter "len" is a length parameter
  int echo_array_with_result(int[] int_in, int[] int_out, int len);
  
  
  // parameter "string_in" is an input parameter
  int print_error_string(String[] string_in);
  
  
  // parameter "string_inout1" is an inout parameter
  // parameter "string_inout2" is an inout parameter
  int echo_strings(String[] string_inout1, String[] string_inout2);
  
  
  // parameter "double_in" is an input parameter
  // parameter "double_out" is an output parameter
  int echo_double(double[] double_in, double[] double_out);
  
  
  // parameter "float_in" is an input parameter
  // parameter "float_out" is an output parameter
  int echo_float(float[] float_in, float[] float_out);
  
  
  // parameter "in_out" is an inout parameter
  // parameter "len" is a length parameter
  int echo_inout_array_with_result(int[] in_out, int len);
  
  
  // parameter "int_in" is an input parameter
  // parameter "int_out" is an output parameter
  // parameter "len" is a length parameter
  void echo_array(int[] int_in, int[] int_out, int len);
  
  
  // parameter "in" is an input parameter
  // parameter "out" is an output parameter
  int echo_long_long_int(long[] in, long[] out);
  
  
  // parameter "string_in" is an input parameter
  // parameter "string_out" is an output parameter
  int echo_string(String[] string_in, String[] string_out);
  
  
  // parameter "input" is an input parameter
  // parameter "output" is an output parameter
  int echo_logical(boolean[] input, boolean[] output);
  
  
  // parameter "string_in" is an input parameter
  int print_string(String[] string_in);
  
  
  // parameter "int_in1" is an input parameter
  // parameter "int_in2" is an input parameter
  // parameter "int_out1" is an output parameter
  // parameter "int_out2" is an output parameter
  // parameter "len" is a length parameter
  int echo_2_int(int[] int_in1, int[] int_in2, int[] int_out1, int[] int_out2, int len);
  
  
  // parameter "int_in" is an input parameter
  // parameter "int_out" is an output parameter
  int echo_int_array(int[] int_in, int[] int_out);
  
}
