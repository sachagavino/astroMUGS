class Star:
    def __init__(self, mass=0.5, luminosity=1.0, temperature=4000., x=0.0, y=0.0, z=0.0):
        self.mass = mass
        self.luminosity = luminosity
        self.temperature = temperature
        self.radius = (5780/temperature)**2*(luminosity)**(1/2)
        self.x = x
        self.y = y
        self.z = z