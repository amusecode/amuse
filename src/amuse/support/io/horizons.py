from datetime import date, timedelta

class LoadStar(object):

    def __init__(self, planet):
        f = open('../viewer/'+planet+'.txt','r')
        M = f.readlines()
        f.close()
        """
        radiusline = M[4]
        radius = float(radiusline.split('=')[1].split('(')[0])

        massline =  M[5]
        Myunit = massline.split('(')[1].split(')')[0].split('^')
        Unitbase = float(Myunit[0])
        Expo = float(Myunit[1].split('kg')[0])
        mass = Unitbase**Expo
        """
        days = []

        #find '$$SOE' this indicates the start of the R,V, T data
        #we do not expect more than one $$OE entry
                
        start_r_v_data_index = [i for i, line in enumerate(M) if '$$SOE' in line][0]

        for i, s in enumerate(M):
            if 'A.D.' in s and i>start_r_v_data_index:
                days.append(i)
        
        r = []
        v = []
        julian_day = []
        
        for i in days:
            
            julian_day.append((float(M[i].split('=')[0])-1721424.5))

            #http://ssd.jpl.nasa.gov/?horizons_doc#time
            rs = M[i+1].split(' ')
            vs = M[i+2].split(' ')
            foor = []
            foov = []
            for j, n in enumerate(rs):
                if not n == '':
                    foor.append(float(n))
            for j, n in enumerate(vs):
                if not n == '':
                    foov.append(float(n))
            
                    
            r.append([foor[0], foor[1], foor[2]])
            v.append([foov[0], foov[1], foov[2]])
            
        self.ordinal = julian_day
        self.r = r
        self.v = v
        #not implemented yet:
        self.mass = 0.0
        self.radius = 0.0
        self.name = planet
        self.max = self.ordinal[-1]
        del M

    def get_vectors_at_date(self, at_date):
        if at_date.toordinal()<self.max:
            indices = [i for i, j in enumerate(self.ordinal) if j==at_date.toordinal()]
            if len(indices)==0:
                print "no ordinal index found"
                return [0,0,0],[0,0,0]
            else:
                return self.r[indices[0]], self.v[indices[0]]
        else:
            print "horizons indexing ERROR"
            return self.r[-1],self.v[-1]

class NewStar(LoadStar):
    def __init__(self, *kargs):
        self.ordinal = []
        self.r = []
        self.v = []
        self.mass = 0.0
        self.radius = 0.0
        self.name = ""
        self.max = 0



    
