package ibis.amuse;

/**
 * Exception cause by some error in either starting or running an amuse code.
 * 
 * @author Niels Drost
 *
 */
public class CodeException extends Exception {

    private static final long serialVersionUID = 1L;

    public CodeException() {
        super();
    }

    public CodeException(String message, Throwable cause) {
        super(message, cause);
    }

    public CodeException(String message) {
        super(message);
    }

    public CodeException(Throwable cause) {
        super(cause);
    }

}
