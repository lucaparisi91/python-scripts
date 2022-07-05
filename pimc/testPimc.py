import unittest
from pimc import *
import json
import numpy as np

class TestParameter(unittest.TestCase):

    def testLoad(self):
        with open("input.json") as f:
            j=json.load(f)
        
        particles = jSonParameter("particles")

        assert(particles[j] == [100] )

        N = nParticlesParameter()

        assert(N[j]==100)

        V=volumeParameter()

        np.testing.assert_almost_equal( 12.915496650148839**3, V[j]  )


        V0 = jSonParameter("action[0]/potential/V0",label="V0")
        np.testing.assert_almost_equal(V0[j] , 11.057432794419006 )

        alpha = jSonParameter("action[0]/potential/alpha",label="alpha")
        np.testing.assert_almost_equal(6.4577483250744185 , alpha[j]   )

        nBeads = jSonParameter("nBeads",label="M")

        self.assertEqual(nBeads[j] , 10 )
        

















if __name__ == '__main__':
    unittest.main()

