from math import log10
from classes import Particle
from vpython import vector, rate, mag, sqrt

dt = 1e-23


def update_forces(particles):
    closest = 1e99
    #en_sum = 0
    for part in particles:
        part.set_force(particles)
        if part.closest < closest:
            closest = part.closest
        #pot = part.get_potential_energy(particles)
        #en_sum += part.get_kinetic_energy() + pot
        # print(part.get_kinetic_energy(), " : ", pot,
        #      " : ", part.get_kinetic_energy() + pot)

    for part in particles:
        if part.collided:
            if part.collision(part.collided):
                particles.append(Particle(identity=len(particles), particle_type="e+", pos=part.pos,
                                 velocity=vector(sqrt(2*part.get_kinetic_energy()/9.10938291e-31), 0, 0)))
                particles[-1].pos.mag += part.radius + particles[-1].radius
                part.velocity_next = vector(0, 0, 0)
    # print(en_sum)
    return closest


def update(particles):
    global dt
    closest = update_forces(particles)
    fastest = 0
    for part in particles:
        part.set_velocity(dt)
        part.set_pos(dt)
        if mag(part.velocity) > fastest:
            fastest = mag(part.velocity)
        # print(particles.index(part), " : ", part.force,
        #      " : ", part.velocity, " : ", part.pos)

    # print(dt)
    if closest > 0:
        dt = 10**(-10.1249999992 * 0.94408751129 **
                  # sqrt(max(log10(fastest), 0) + 5) * 0.28571428571)
                  (log10(closest)) * (0.09375*max(log10(fastest), 0)+.25))


def main():
    particles = []
    for i in range(2):
        particles.append(Particle(identity=i, particle_type="p", pos=vector(
            i*1e-14 - 1e-16 * i, 0, 0), velocity=vector((0 if i == 1 else 1)*100000000 * (-1)**int(.5*i), 0, 0)))
    # particles.append(electron(identity=2, pos=vector(
    #    -5e-14, 1e-12, 0), velocity=vector(500000, 50000, 200000)))

    while True:
        rate(500)
        # input()
        update(particles)


if __name__ == "__main__":
    main()
