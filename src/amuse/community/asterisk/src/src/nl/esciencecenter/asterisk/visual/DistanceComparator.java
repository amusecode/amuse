package nl.esciencecenter.asterisk.visual;

import java.util.Comparator;

import nl.esciencecenter.esight.math.VecF3;
import nl.esciencecenter.esight.math.VectorFMath;

public class DistanceComparator implements Comparator<StarModel> {
    VecF3 cameraPosition;

    public DistanceComparator(VecF3 cameraPosition) {
        this.cameraPosition = cameraPosition;
    }

    @Override
    public int compare(StarModel m1, StarModel m2) {
        VecF3 m1Location = m1.getLocation();
        VecF3 m2Location = m2.getLocation();

        if (VectorFMath.length(m1Location.sub(cameraPosition)) > VectorFMath
                .length(m2Location.sub(cameraPosition))) {
            return -1;
        } else if (VectorFMath.length(m1Location.sub(cameraPosition)) < VectorFMath
                .length(m2Location.sub(cameraPosition))) {
            return 1;
        } else {
            return 0;
        }
    }
}
