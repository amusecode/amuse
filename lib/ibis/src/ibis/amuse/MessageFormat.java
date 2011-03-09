package ibis.amuse;

/**
 * Some constants describing the message format. A Message (either a call or a
 * result) consists of a header describing the call id, function id, number of
 * calls in this message, and number of parameters of each type. Types are send
 * as their byte representation, except boolean which is send as a single byte,
 * and strings which are send as UTF-8. Also, the strings are preceded by an
 * additional header listing the size of all following strings.
 * 
 * @author Niels Drost
 * 
 */
public class MessageFormat {

    public static final String MAGIC_STRING = "magic_string";
    
    public static final int HEADER_SIZE = 9; // integers

    public static final int HEADER_CALL_ID_INDEX = 0;
    public static final int HEADER_FUNCTION_ID_INDEX = 1;
    public static final int HEADER_CALL_COUNT_INDEX = 2;
    public static final int HEADER_INT_COUNT_INDEX = 3;
    public static final int HEADER_LONG_COUNT_INDEX = 4;
    public static final int HEADER_FLOAT_COUNT_INDEX = 5;
    public static final int HEADER_DOUBLE_COUNT_INDEX = 6;
    public static final int HEADER_BOOLEAN_COUNT_INDEX = 7;
    public static final int HEADER_STRING_COUNT_INDEX = 8;

    public static final int SIZEOF_INT = 4;
    public static final int SIZEOF_LONG = 8;
    public static final int SIZEOF_FLOAT = 8;
    public static final int SIZEOF_DOUBLE = 4;
    public static final int SIZEOF_BOOLEAN = 1;

    public static final byte TRUE = (1 & 0xFF);
    public static final byte FALSE = (0 & 0xFF);

    public static final int FUNCTION_ID_INIT = 10101010;
    public static final int FUNCTION_ID_STOP = 0;
    
    public static final int FUNCTION_ID_REDIRECT_OUTPUT = 1141573512;

}
