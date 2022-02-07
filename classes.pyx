import cython
from math import log10
from vpython import sphere, arrow, label, vector, hypot, color, mag2

class Particle(sphere):
    def __init__(self, identity, particle_type, pos, velocity = vector(0,0,0)):
        self.identity = identity
        self.particle_type = particle_type
        self.labl = label(pos=pos, text=f"{particle_type}:{self.identity}")
        self.force = vector(0, 0, 0) # Kgms^-2
        self.velocity = velocity # ms^-1
        self.velocity_next = None
        self.closest = 1e99
        self.collided = None

        if particle_type == "p":
            super().__init__(pos=pos, radius=8.4e-16, color=color.red)
            self.mass = 1.672621777e-27  # kg
            self.charge = 1.602176565e-19  # C
        elif particle_type == "n":
            super().__init__(pos=pos, radius=8.4e-16, color=color.white)
            self.mass = 1.67493e-27  # kg
            self.charge = 0  # C
        elif particle_type == "e-":
            super().__init__(pos=pos, radius=2.8e-16, color=color.blue)
            self.mass = 9.10938291e-31  # kg
            self.charge = -1.602176565e-19  # C
        elif particle_type == "e+":
            super().__init__(pos=pos, radius=2.8e-16, color=color.yellow)
            self.mass = 9.10938291e-31  # kg
            self.charge = 1.602176565e-19  # C
    
    def collision(self, part):
        self.collided = None
        # Solving for v in completely elastic equations
        vel_after_collision_self = (self.velocity+part.mass/self.mass*(2*part.velocity - self.velocity))/(1+part.mass/self.mass)
        vel_after_collision_part = self.mass/part.mass*(self.velocity - vel_after_collision_self) + part.velocity
        self.velocity_next = vel_after_collision_self
        print(f"boom {self.identity} (self after, sum after, other after collision):", vel_after_collision_self, self.velocity_next, vel_after_collision_part)
        if self.particle_type == "p" and part.particle_type == "p" and self.identity > part.identity:
            self.particle_type = "n"
            self.labl.text = f"n:{self.identity}"
            self.mass = 1.67493e-27  # kg
            self.charge = 0  # C
            self.color = color.white
            return True
        return False

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
                    if self.collided:
                        print("Simultatious collision! WARNING: This is not a correct representation of what is supposed to happen!")
                    self.collided = part

                one_over_r = 1 / r
                one_over_r2 = one_over_r**2
                f_el = f*self.charge * part.charge * one_over_r2
                self.force += r_vec * f_el * one_over_r * -1
                if r < 1e-14:
                    f_g = G*self.mass * part.mass * one_over_r2
                    #print("strong")
                    # Strong force crude aproximation attempt (wikipedia: At the range of 10^−15 m,
                    # the strong force is approximately 10^38 times as strong as gravitation.
                    f_s = 10**38*f_g#(10**38/(r*10**15))*f_g
                    #f_s = 7.6388794e-98 * 2.63157895e-7 ** (log10(r)) * f_g
                    self.force += r_vec * f_s * one_over_r
                
        return self.force

    def get_kinetic_energy(self):
        #print(f".5 * self.mass * mag2(self.velocity) = .5 * {self.mass} * {mag2(self.velocity)} = .5 * {self.mass} * {self.velocity.x**2}")
        return .5 * self.mass * mag2(self.velocity)

    def get_potential_energy(self, particles):
        G = 6.67384e-11  # Nm^2kg^-2
        f = 8.987551787e9  # Nm^2C^-2
        potential_energy = 0
        for part in particles:
            if part.identity != self.identity:
                r = hypot(part.pos.x-self.pos.x, part.pos.y-self.pos.y,
                        part.pos.z-self.pos.z)
                one_over_r = 1 / r
                e_el = f*self.charge * part.charge * one_over_r
                print(f"r: {r}, f*self.charge * part.charge * one_over_r = f*{self.charge * part.charge} * {one_over_r} = {e_el}")
                potential_energy += e_el
                if r < 1e-14:
                    e_g = -G*self.mass * part.mass * one_over_r
                    #print("strong")
                    # Strong force crude aproximation attempt (wikipedia: At the range of 10^−15 m,
                    # the strong force is approximately 10^38 times as strong as gravitation.
                    e_s = 10**38*e_g#(10**38/(r*10**15))*e_g
                    #f_s = 7.6388794e-98 * 2.63157895e-7 ** (log10(r)) * f_g
                    potential_energy += e_s
                
        return potential_energy

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

