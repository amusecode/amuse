import nl.esciencecenter.asterisk.*;

public class CodeParticle {

    public static final int TYPE_STAR = 0;
    public static final int TYPE_POINT_GAS = 1;
    public static final int TYPE_OCTREE_GAS = 2;
    public static final int TYPE_SPHERE = 3;
    public static final int TYPE_MARKER = 4;

    private int index;
    private double x;
    private double y;
    private double z;
    private double radius;
    private double red;
    private double green;
    private double blue;
    private double opacity;
    private int type;

    CodeParticle(int index,  double x, double y, double z, double radius, double red, double green,
            double blue, double opacity, int type) {
        this.index = index;
        this.x = x;
        this.y = y;
        this.z = z;
        this.radius = radius;
        this.red = red;
        this.green = green;
        this.blue = blue;
        this.opacity = opacity;
        this.type = type;
    }

    PointGas asPointGas() {
        float[] coordinates = new float[] { (float) x, (float) y, (float) z };

        float[] color = new float[] { (float) red, (float) green, (float) blue, (float) opacity };
        
        return new PointGas(index, coordinates, color);
    }

    Sphere asSphere() {
        float[] coordinates = new float[] { (float) x, (float) y, (float) z };

        float[] color = new float[] { (float) red, (float) green, (float) blue, (float) opacity };

        return new Sphere(index, coordinates, (float) radius, color);
    }

    SPHGas asSphGas() {
        float[] coordinates = new float[] { (float) x, (float) y, (float) z };

        float[] color = new float[] { (float) red, (float) green, (float) blue, (float) opacity };

        return new SPHGas(index, coordinates, color);
    }

    Star asStar() {
        float[] coordinates = new float[] { (float) x, (float) y, (float) z };

        float[] color = new float[] { (float) red, (float) green, (float) blue, (float) opacity };
        
        return new Star(index, coordinates, (float) radius, color);
    }

    public double getOpacity() {
        return opacity;
    }

    public double getBlue() {
        return blue;
    }

    public double getGreen() {
        return green;
    }

    public int getIndex() {
        return index;
    }

    public double getRadius() {
        return radius;
    }

    public double getRed() {
        return red;
    }

    public int getType() {
        return type;
    }

    public double getX() {
        return x;
    }

    public double getY() {
        return y;
    }

    public double getZ() {
        return z;
    }

    public void setOpacity(double opacity) {
        this.opacity = opacity;
    }

    public void setBlue(double blue) {
        this.blue = blue;
    }

    public void setColor(double red, double green, double blue) {
        this.red = red;
        this.green = green;
        this.blue = blue;
    }

    public void setGreen(double green) {
        this.green = green;
    }

    public void setIndex(int index) {
        this.index = index;
    }

    public void setPosition(double x, double y, double z) {
        this.x = x;
        this.y = y;
        this.z = z;
        
    }

    public void setRadius(double radius) {
        this.radius = radius;
    }

    public void setRed(double red) {
        this.red = red;
    }

    public void setType(int type) {
        this.type = type;
        
    }

    public void setX(double x) {
        this.x = x;
    }

    public void setY(double y) {
        this.y = y;
    }

    public void setZ(double z) {
        this.z = z;
    }

    @Override
    public String toString() {
        return "CodeParticle [index=" + index + ", type=" + type + ", x=" + x + ", y=" + y + ", z=" + z + ", radius="
                + radius + ", red=" + red + ", green=" + green + ", blue=" + blue + ", opacity=" + opacity + "]";
    }
    
    
}
