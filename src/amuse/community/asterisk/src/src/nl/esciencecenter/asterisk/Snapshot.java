package nl.esciencecenter.asterisk;

import java.io.Serializable;

public class Snapshot implements Serializable {
    private static final long serialVersionUID = 4342001313150488429L;

    private final Sphere[] spheres;

    private final Star[] stars;

    private final SPHGas[] sphGas;
    private final PointGas[] pointGas;

    private final String decription;

    public Snapshot(String decription, Sphere[] spheres, Star[] stars, SPHGas[] sphGas, PointGas[] pointGas) {
        this.spheres = spheres;
        this.stars = stars;
        this.sphGas = sphGas;
        this.pointGas = pointGas;
        this.decription = decription;
    }

    public Sphere[] getSpheres() {
        return spheres;
    }

    public Star[] getStars() {
        return stars;
    }

    public SPHGas[] getSphGas() {
        return sphGas;
    }

    public PointGas[] getPointGas() {
        return pointGas;
    }

    public String getSceneDecriptionString() {
        return decription;
    }
}
