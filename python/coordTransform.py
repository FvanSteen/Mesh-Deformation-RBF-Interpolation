import numpy as np

class coordTransform:
    def toCylindrical(v):
        v_new = np.empty(np.shape(v), dtype=np.double)
        v_new[:,0] = np.sqrt( v[:,0]**2 + v[:,1]**2)
        v_new[:,1] = np.arctan2(v[:,1],v[:,0])
        v_new[:,2] = v[:,2]
        return v_new