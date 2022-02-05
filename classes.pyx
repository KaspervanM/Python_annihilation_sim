import cython
from vpython import sphere, label, vector, hypot


class proton(sphere):
    def __init__(self, identity, pos, velocity = vector(0,0,0)):
        super().__init__(pos=pos, radius=8.4e-16)
        self.identity = identity
        self.force = vector(0, 0, 0) # Kgms^-2
        self.velocity = velocity # ms^-1
        self.mass = 1.672621777e-27  # kg
        self.charge = 1.602176565e-19  # C
        self.labl = label(pos=self.pos, text=f"p:{self.identity}")

    def set_force(self, particles):
        G = 6.67384e-11  # Nm^2kg^-2
        f = 8.987551787e9  # Nm^2C^-2
        self.force = vector(0,0,0)
        for part in particles:
            if part.identity != self.identity:
                r_vec = vector(self.pos.x-part.pos.x, self.pos.y -
                           part.pos.y, self.pos.z-part.pos.z)
                one_over_r = 1 / hypot(self.pos.x-part.pos.x, self.pos.y -
                        part.pos.y, self.pos.z-part.pos.z)
                one_over_r2 = one_over_r**2
                f_g = G*self.mass * part.mass * one_over_r2
                f_em = f*self.charge * part.charge * one_over_r2
                self.force += r_vec * f_g * one_over_r
                self.force += r_vec * f_em * one_over_r
        return self.force

    def get_momentum(self):
        return self.mass*self.velocity

    def get_acceleration(self):
        return self.force / self.mass

    def set_velocity(self):
        self.velocity += self.get_acceleration()
        return self.velocity

    def set_pos(self, dt):
        self.pos += self.velocity * dt
        self.labl.pos = self.pos
        return self.pos