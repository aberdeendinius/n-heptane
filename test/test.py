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

species(
    label='Ar',
    reactive=False,
    structure=SMILES("[Ar]"),
)


simpleReactor(
    temperature=(1600,'K'),
    pressure=(400,'Pa'),
    initialMoleFractions={
        "n-heptane": 0.02,
        "Ar": 0.98,
    },
    terminationConversion={
        'n-heptane': 0.99,
    },
    terminationTime=(1e6,'s'),
)

simpleReactor(
    temperature=(2000,'K'),
    pressure=(400,'Pa'),
    initialMoleFractions={
        "n-heptane": 0.02,
        "Ar": 0.98,
    },
    terminationConversion={
        'n-heptane': 0.99,
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
    temperatures=(300,3000,'K',8),
    pressures=(0.001,100,'bar',5),
    interpolation=('Chebyshev', 6, 4),
)
