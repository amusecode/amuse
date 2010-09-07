FUNCTION distance_between_points(x0,y0,z0,x1,y1,z1)
    DOUBLE PRECISION :: distance_between_points
    DOUBLE PRECISION :: x0,y0,z0,x1,y1,z1
    DOUBLE PRECISION :: dx,dy,dz
    dx = x1 - x0;
    dy = y1 - y0;
    dz = z1 - z0;
    distance_between_points = sqrt(dx * dx  + dy * dy + dz * dz);
END FUNCTION


FUNCTION find_nearest_neighbors(npoints,x,y,z,n0,n1,n2)
    INTEGER :: npoints
    INTEGER :: i,j, k, kk
    INTEGER :: find_nearest_neighbors
    DOUBLE PRECISION :: x(npoints),y(npoints),z(npoints)
    INTEGER :: n0(npoints),n1(npoints),n2(npoints)
    DOUBLE PRECISION :: r
    INTEGER :: nn_index(3)
    DOUBLE PRECISION :: nn_distance(3)
    
    DO i = 1, npoints, 1
        DO k = 1, 3
            nn_index(k) = 0
            nn_distance(k) = 0.0
        END DO
        
        DO j = 1, npoints, 1
            IF (i.NE.j) THEN
                r = distance_between_points(x(i), y(i), z(i), &
                    x(j), y(j), z(j))
                
                
                DO k = 1, 3
                    IF (nn_index(k).EQ.0) THEN
                        nn_index(k) = j
                        nn_distance(k) = r
                    ELSE
                        IF (r.LT.nn_distance(k)) THEN
                            DO kk = 3, k, -1
                                nn_index(kk) = nn_index(kk-1)
                                nn_distance(kk) = nn_distance(kk-1)
                            END DO
                            nn_index(k) = k
                            nn_distance(k) = r
                            EXIT
                        END IF
                    END IF
                    nn_distance(k) = 0.0
                END DO
            END IF
        END DO
        
        n0(i) = nn_index(1)
        n1(i) = nn_index(2)
        n2(i) = nn_index(3)
    END DO
    
    find_nearest_neighbors = 0
END FUNCTION
