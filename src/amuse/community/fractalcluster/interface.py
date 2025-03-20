from amuse.support.import_helper import load_code


FractalClusterInterface = load_code("fractalcluster", "FractalClusterInterface")
FractalCluster = load_code("fractalcluster", "FractalCluster")
MakeFractalCluster = load_code("fractalcluster", "MakeFractalCluster")
new_fractal_cluster_model = load_code("fractalcluster", "new_fractal_cluster_model")

Fractalcluster = FractalCluster
