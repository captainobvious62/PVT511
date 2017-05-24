import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt

# Increase the starting pressure a bit, behavior at very low pressure is problematic
CP.set_config_double(CP.PHASE_ENVELOPE_STARTING_PRESSURE_PA, 1e4)

SRK = CP.AbstractState('SRK','Methane&Ethane')
SRK.set_mole_fractions([0.5, 1 - 0.5])
for kij, c in zip([0.0, 0.1],['r','b']):

    # Set the interaction parameter
    SRK.set_binary_interaction_double(0, 1, "kij", kij)

    # Some VLE calculations
    for p in [1e4, 1e5, 1e6]:
        SRK.update(CP.PQ_INPUTS, p, 0)
        plt.plot(SRK.T(), SRK.p(), '<', color = c)

        SRK.update(CP.PQ_INPUTS, p, 1)
        plt.plot(SRK.T(), SRK.p(), '>', color = c)

    # Phase envelope
    SRK.build_phase_envelope("")
    PE = SRK.get_phase_envelope_data()
    plt.plot(PE.T, PE.p, '-', label = '$k_{ij} = $' + str(kij), color = c)

    # Critical point
    pts = SRK.all_critical_points()
    for pt in pts:
      plt.plot(pt.T, pt.p, '*', color = c)

# A phase envelope calculated with SRK transformations in a multi-fluid model
HEOS = CP.AbstractState('HEOS','Methane-SRK&Ethane-SRK')
HEOS.set_mole_fractions([0.5, 0.5])
HEOS.build_phase_envelope("none")
PE = HEOS.get_phase_envelope_data()
plt.plot(PE.T, PE.p, '-', label = 'SRK with transformations in multi-fluid', color = 'g')

plt.xlabel('Temperature [K]')
plt.ylabel('Pressure [Pa]')
plt.yscale('log')
plt.legend(loc='best')
plt.tight_layout()
plt.show()
