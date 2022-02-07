from math import log10
from classes import Proton, Electron
from vpython import vector, rate, mag, sqrt

dt = 1e-23


def update_forces(particles):
    closest = 1e99
    for part in particles:
        part.set_force(particles)
        if part.closest < closest:
            closest = part.closest
    return closest


def update(particles):
    global dt
    closest = update_forces(particles)
    fastest = 0
    for part in particles:
        part.set_velocity(dt)
        if mag(part.velocity) > fastest:
            fastest = mag(part.velocity)
        part.set_pos(dt)
        # print(particles.index(part), " : ", part.force,
        #      " : ", part.velocity, " : ", part.pos)

    #   print(dt, closest, fastest)
    if closest > 0 and fastest > 1:
        dt = 10**(-7.15743510167 * 1.08695652174 **
                  (-log10(closest)) * sqrt(log10(fastest))/2.72)


def main():
    particles = []
    for i in range(2):
        particles.append(Proton(identity=i, pos=vector(
            i*5e-14 - 1e-16 * i, 0, 0), velocity=vector((0 if i == 1 else 1)*100000000 * (-1)**int(.5*i), 0, 0)))
    # particles.append(electron(identity=2, pos=vector(
    #    -5e-14, 1e-12, 0), velocity=vector(500000, 50000, 200000)))

    while True:
        rate(300)
        # input()
        update(particles)


if __name__ == "__main__":
    main()
