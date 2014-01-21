package nl.esciencecenter.asterisk.input;

import nl.esciencecenter.esight.input.InputHandler;
import nl.esciencecenter.esight.math.VecF3;

import com.jogamp.newt.event.KeyEvent;
import com.jogamp.newt.event.KeyListener;
import com.jogamp.newt.event.MouseEvent;
import com.jogamp.newt.event.MouseListener;

/**
 * @author Maarten van Meersbergen <m.van.meersbergen@esciencecenter.nl> A
 *         singleton pattern generic Input event Handler for use in OpenGL
 *         applications. Currently handles only basic mouse events
 *         (left-click-drag, scrollwheel).
 * 
 */
public class AsteriskInputHandler extends InputHandler implements MouseListener, KeyListener {

    /** Initial value for the rotation in the X direction */
    protected float rotationXorigin = 0;
    /**
     * Final rotation in the X direction, translated to openGL units, stored to
     * make successive rotations smooth
     */
    protected float rotationX;

    /** Initial value for the rotation in the Y direction */
    protected float rotationYorigin = 0;
    /**
     * Final rotation in the Y direction, translated to openGL units, stored to
     * make successive rotations smooth
     */
    protected float rotationY;

    /** Initial value for the rotation in the Y direction */
    protected float rotationZorigin = 0;
    /**
     * Final rotation in the Y direction, translated to openGL units, stored to
     * make successive rotations smooth
     */
    protected float rotationZ;

    /** Mouse drag start point in X direction */
    protected float dragLeftXorigin;
    /** Mouse drag start point in Y direction */
    protected float dragLeftYorigin;

    /** Final rotation in openGL units */
    public VecF3 rotation, translation;
    /** Final view distance (translation) in openGL units */
    public float viewDist = -3f;
    /** Current direction of the view */
    private octants current_view_octant = octants.PPP;

    private float translationX = 0f;
    private float translationY = 0f;

    private float translationXorigin = 0f;
    private float translationYorigin = 0f;

    private static class SingletonHolder {
        public static final AsteriskInputHandler instance = new AsteriskInputHandler();
    }

    /**
     * The only access point for this singleton class.
     * 
     * @return The only instance of this class allowed at one time.
     */
    public static AsteriskInputHandler getInstance() {
        return SingletonHolder.instance;
    }

    protected AsteriskInputHandler() {
        reset();
    }

    private void reset() {
        rotationXorigin = 0;
        rotationX = 0;
        rotationYorigin = 0;
        rotationY = 0;
        rotationZorigin = 0;
        rotationZ = 0;

        dragLeftXorigin = 0;
        dragLeftYorigin = 0;

        translationX = 0f;
        translationXorigin = 0f;
        translationY = 0f;
        translationYorigin = 0f;

        rotation = new VecF3();
        translation = new VecF3();
        viewDist = -3f;
        current_view_octant = octants.PPP;
    }

    @Override
    public void mouseClicked(MouseEvent e) {
    }

    @Override
    public void mouseEntered(MouseEvent e) {
        // Empty - unneeded
    }

    @Override
    public void mouseExited(MouseEvent e) {
        // Empty - unneeded
    }

    @Override
    public void mousePressed(MouseEvent e) {
        // if (e.isButtonDown(MouseEvent.BUTTON1) ||
        // e.isButtonDown(MouseEvent.BUTTON3)) {
        dragLeftXorigin = e.getX();
        dragLeftYorigin = e.getY();
        // }
    }

    @Override
    public void mouseReleased(MouseEvent e) {
        rotationXorigin = rotationX;
        rotationYorigin = rotationY;
        rotationZorigin = rotationZ;

        translationXorigin = translationX;
        translationYorigin = translationY;
    }

    @Override
    public void mouseDragged(MouseEvent e) {
        if (e.isButtonDown(MouseEvent.BUTTON1)) {
            // x/y reversed because of axis orientation. (up/down => x axis
            // rotation in OpenGL)
            if (e.isControlDown()) {
                rotationZ = ((e.getX() - dragLeftXorigin) + rotationZorigin) % 360;
            } else {
                if (e.isShiftDown()) {
                    rotationX = ((e.getX() - dragLeftXorigin) / 10f + rotationXorigin) % 360;
                    rotationY = ((e.getY() - dragLeftYorigin) / 10f + rotationYorigin) % 360;
                } else {
                    rotationX = ((e.getX() - dragLeftXorigin) + rotationXorigin) % 360;
                    rotationY = ((e.getY() - dragLeftYorigin) + rotationYorigin) % 360;
                }
            }
            // Make sure the numbers are always positive (so we can determine
            // the octant we're in more easily)
            if (rotationX < 0) {
                rotationX = 360f + rotationX % 360;
            }
            if (rotationY < 0) {
                rotationY = 360f + rotationY % 360;
            }
            if (rotationZ < 0) {
                rotationZ = 360f + rotationZ % 360;
            }

            rotation.set(0, rotationY);
            rotation.set(1, rotationX);
            rotation.set(2, rotationZ);
            setCurrentOctant(rotation);
        } else if (e.isButtonDown(MouseEvent.BUTTON3)) {
            if (e.isShiftDown()) {
                translationX = (.01f * (e.getX() - dragLeftXorigin)) / 10f + translationXorigin;
                translationY = (-.01f * (e.getY() - dragLeftYorigin)) / 10f + translationYorigin;
            } else {
                translationX = (.01f * (e.getX() - dragLeftXorigin)) + translationXorigin;
                translationY = (-.01f * (e.getY() - dragLeftYorigin)) + translationYorigin;
            }
            translation.set(0, translationX);
            translation.set(1, translationY);
        }
    }

    @Override
    public void mouseMoved(MouseEvent e) {
        // Empty - unneeded
    }

    @Override
    public void mouseWheelMoved(MouseEvent e) {
        float newViewDist = this.viewDist;

        if (e.isShiftDown()) {
            newViewDist += e.getWheelRotation() * 0.02f;
        } else {
            newViewDist += e.getWheelRotation() * 0.10f;
        }
        viewDist = newViewDist;
    }

    /**
     * Setter for the current viewing 'angle' by octant.
     * 
     * @param rotation
     *            The rotation from which to calculate the viewing octant.
     */
    private void setCurrentOctant(VecF3 rotation) {
        float x = rotation.get(0);
        int qx = (int) Math.floor(x / 90f);
        float y = rotation.get(1);
        int qy = (int) Math.floor(y / 90f);

        if (qx == 0 && qy == 0) {
            current_view_octant = octants.NPP;
        } else if (qx == 0 && qy == 1) {
            current_view_octant = octants.NPN;
        } else if (qx == 0 && qy == 2) {
            current_view_octant = octants.PPN;
        } else if (qx == 0 && qy == 3) {
            current_view_octant = octants.PPP;

        } else if (qx == 1 && qy == 0) {
            current_view_octant = octants.PPN;
        } else if (qx == 1 && qy == 1) {
            current_view_octant = octants.PPP;
        } else if (qx == 1 && qy == 2) {
            current_view_octant = octants.NPP;
        } else if (qx == 1 && qy == 3) {
            current_view_octant = octants.NPN;

        } else if (qx == 2 && qy == 0) {
            current_view_octant = octants.PNN;
        } else if (qx == 2 && qy == 1) {
            current_view_octant = octants.PNP;
        } else if (qx == 2 && qy == 2) {
            current_view_octant = octants.NNP;
        } else if (qx == 2 && qy == 3) {
            current_view_octant = octants.NNN;

        } else if (qx == 3 && qy == 0) {
            current_view_octant = octants.NNP;
        } else if (qx == 3 && qy == 1) {
            current_view_octant = octants.NNN;
        } else if (qx == 3 && qy == 2) {
            current_view_octant = octants.PNN;
        } else if (qx == 3 && qy == 3) {
            current_view_octant = octants.PNP;
        }
    }

    @Override
    public void keyPressed(KeyEvent arg0) {
        if (arg0.getKeyCode() == KeyEvent.VK_SPACE) {
            reset();
        }

    }

    @Override
    public void keyReleased(KeyEvent arg0) {
        if (arg0.getKeyCode() == KeyEvent.VK_SPACE) {
            reset();
        }
    }

    public void keyTyped(KeyEvent arg0) {
        if (arg0.getKeyCode() == KeyEvent.VK_SPACE) {
            reset();
        }
    }

    /**
     * 
     * @return the current OpenGL ModelView rotation variable
     */
    @Override
    public VecF3 getRotation() {
        return rotation;
    }

    /**
     * @param rotation
     *            the current OpenGL ModelView rotation variable to set
     */
    @Override
    public void setRotation(VecF3 rotation) {
        // reset();
        this.rotation = rotation;
    }

    /**
     * 
     * @return the current OpenGL ModelView translation variable
     */
    public VecF3 getTranslation() {
        return translation;
    }

    /**
     * @param rotation
     *            the OpenGL ModelView translation variable to set
     */
    public void setTranslation(VecF3 translation) {
        this.translation = translation;
    }

    /**
     * @return the current OpenGL ModelView View distance (translation) variable
     */
    @Override
    public float getViewDist() {
        return viewDist;
    }

    /**
     * @param viewDist
     *            the current OpenGL ModelView View distance (translation)
     *            variable to set
     */
    @Override
    public void setViewDist(float viewDist) {
        this.viewDist = viewDist;
    }

    /**
     * @return the current_view_octant (mainly used for octrees)
     */
    @Override
    public octants getCurrentViewOctant() {
        return current_view_octant;
    }
}
