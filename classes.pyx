import cython
from math import log10
from vpython import sphere, arrow, label, vector, hypot, color


class proton(sphere):
    def __init__(self, identity, pos, velocity = vector(0,0,0)):
        super().__init__(pos=pos, radius=8.4e-16, color=color.red)
        self.identity = identity
        self.force = vector(0, 0, 0) # Kgms^-2
        self.velocity = velocity # ms^-1
        self.velocity_next = None
        self.mass = 1.672621777e-27  # kg
        self.charge = 1.602176565e-19  # C
        self.labl = label(pos=self.pos, text=f"p:{self.identity}")
        self.closest = 1e99

    def set_force(self, particles):
        G = 6.67384e-11  # Nm^2kg^-2
        f = 8.987551787e9  # Nm^2C^-2
        self.force = vector(0,0,0)
        for part in particles:
            if part.identity != self.identity:
                r_vec = vector(part.pos.x-self.pos.x,
                            part.pos.y-self.pos.y, part.pos.z-self.pos.z)
                r = hypot(part.pos.x-self.pos.x, part.pos.y-self.pos.y,
                        part.pos.z-self.pos.z)
                #print("r", r)
                if r < self.closest:
                    self.closest = r
                if r < self.radius + part.radius:
                    # Solving for v in completely elastic equations
                    vel_after_collision_self = (self.velocity+part.mass/self.mass*(2*part.velocity - self.velocity))/(1+part.mass/self.mass)
                    vel_after_collision_part = self.mass/part.mass*(self.velocity - vel_after_collision_self) + part.velocity
                    if self.velocity_next:
                        print("Simultatious collision! WARNING: This is not a correct representation of what is supposed to happen!")
                        #self.velocity_next += vel_after_collision_self
                    else:
                        self.velocity_next = vel_after_collision_self
                    print(f"boom {self.identity} (self after, sum after, other after collision):", vel_after_collision_self, self.velocity_next, vel_after_collision_part)

                one_over_r = 1 / r
                one_over_r2 = one_over_r**2
                f_g = G*self.mass * part.mass * one_over_r2
                f_em = f*self.charge * part.charge * one_over_r2
                self.force += r_vec * f_g * one_over_r
                self.force += r_vec * f_em * one_over_r * -1
                if r < 1e-14:
                    #print("strong")
                    # Strong force crude aproximation attempt (wikipedia: At the range of 10^−15 m,
                    # the strong force is approximately 10^38 times as strong as gravitation.
                    f_s = (10**38/(r*10**15))*f_g
                    #f_s = 7.6388794e-98 * 2.63157895e-7 ** (log10(r)) * f_g
                    self.force += r_vec * f_s * one_over_r
                
        return self.force

    def get_momentum(self):
        return self.mass*self.velocity

    def get_acceleration(self):
        return self.force / self.mass

    def set_velocity(self, dt):
        if not self.velocity_next:
            self.velocity += self.get_acceleration() * dt
            return self.velocity
        self.velocity = self.velocity_next
        self.velocity_next = None
        return self.velocity

    def set_pos(self, dt):
        self.pos += self.velocity * dt
        self.labl.pos = self.pos
        return self.pos


"""class electron(sphere):
    def __init__(self, identity, pos, velocity = vector(0,0,0)):
        super().__init__(pos=pos, radius=2.8e-16, color=color.blue)
        self.identity = identity
        self.force = vector(0, 0, 0) # Kgms^-2
        self.velocity = velocity # ms^-1
        self.mass = 9.10938291e-31  # kg
        self.charge = -1.602176565e-19  # C
        self.labl = label(pos=self.pos, text=f"e:{self.identity}")
        self.closest = 1e99

    def set_force(self, particles):
        G = 6.67384e-11  # Nm^2kg^-2
        f = 8.987551787e9  # Nm^2C^-2
        self.force = vector(0,0,0)
        self.closest = 1e99
        for part in particles:
            if part.identity != self.identity:
                r_vec = vector(part.pos.x-self.pos.x,
                            part.pos.y-self.pos.y, part.pos.z-self.pos.z)
                r = hypot(part.pos.x-self.pos.x, part.pos.y-self.pos.y,
                        part.pos.z-self.pos.z)
                #print("r", r)
                if r < self.closest:
                    self.closest = r
                if r < self.radius + part.radius:
                    vel_after_collision_self = (self.velocity+part.mass/self.mass*(2*part.velocity - self.velocity))/(1+part.mass/self.mass)
                    #vel_after_collision_part = self.mass/part.mass*(self.velocity - vel_after_collision_self) + part.velocity
                    self.velocity = vel_after_collision_self
                    #print("boom", self.velocity)

                one_over_r = 1 / r
                one_over_r2 = one_over_r**2
                f_g = G*self.mass * part.mass * one_over_r2
                f_em = f*self.charge * part.charge * one_over_r2
                self.force += r_vec * f_g * one_over_r
                self.force += r_vec * f_em * one_over_r * -1
                if r < 1e-14:
                    print("strong")
                    # Strong force crude aproximation attempt (wikipedia: At the range of 10^−15 m,
                    # the strong force is approximately 10^38 times as strong as gravitation.
                    f_s = (10**38/(r*10**15))*f_g
                    self.force += r_vec * f_s * one_over_r
                
        return self.force

    def get_momentum(self):
        return self.mass*self.velocity

    def get_acceleration(self):
        return self.force / self.mass

    def set_velocity(self, dt):
        self.velocity += self.get_acceleration() * dt
        return self.velocity

    def set_pos(self, dt):
        self.pos += self.velocity * dt
        self.labl.pos = self.pos
        return self.pos"""