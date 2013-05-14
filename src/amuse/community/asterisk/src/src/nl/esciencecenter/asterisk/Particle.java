package nl.esciencecenter.asterisk;

import java.io.Serializable;

public abstract class Particle implements Serializable {
    private static final long serialVersionUID = -5115248697850683000L;

    private final int         index;
    private final float[]     coordinates;
    private final float       radius;

    private final float[]     color;

    public Particle(int index, float[] coordinates, float radius, float[] color) {
        this.index = index;
        this.coordinates = coordinates;
        this.radius = radius;
        this.color = color;
    }

    public int getIndex() {
        return index;
    }

    public float[] getCoordinates() {
        return coordinates;
    }

    public float getRadius() {
        return radius;
    }

    public float[] getColor() {
        return color;
    }

}
