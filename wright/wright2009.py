#!bin/python
#wright auto-ignition of n-heptane in closed combustion chamber
# Data sources
database(
    thermoLibraries = ['primaryThermoLibrary'],
    reactionLibraries = [],
    seedMechanisms = [],
    kineticsDepositories = ['training'], 
    kineticsFamilies = 'default',
    kineticsEstimator = 'rate rules',
)

# Constraints on generated species
generatedSpeciesConstraints(
    maximumCarbonAtoms = 7,
)

# List of species
species(
    label='n-heptane',
    structure=SMILES("CCCCCCC"),
)

simpleReactor(
    pressure=(80,'bar'),
    temperature=(673, 'K'),
    initialMoleFractions={
        "n-heptane": 1.0,
    },
    terminationTime=(1e6,'s'),
)

#absolute and relative tolerance for ODE solver, atol usually 1e-15 - 1e-25, rtol 1e-4 - 1e-8
simulator(
    atol=1e-16,
    rtol=1e-8,
)

#tolerancemovetocore is how high edge flux ratio of species must reach before entering core model
#toleranceInterruptSimulation how high edge flux ratio to interrupt sim, set equal to tolerancemove to core
model(
    toleranceMoveToCore=0.1,
    toleranceInterruptSimulation=0.1
)

# generates k(T, P) interpolation model
#wo methods available: 'modified strong collision' is faster and less accurate than 'reservoir state'
pressureDependence(
    method='modified strong collision',
    maximumGrainSize=(0.5,'kcal/mol'),
    minimumNumberOfGrains=250,
    temperatures=(650,800,'K',10),
    pressures=(70,250,'bar',10),
    interpolation=('Chebyshev', 6, 4),
)
