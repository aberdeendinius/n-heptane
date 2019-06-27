species(
    label = '[CH]=[C]C=C([CH2])[CH2](19662)',
    structure = SMILES('[CH]=[C]C=C([CH2])[CH2]'),
    E0 = (701.044,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,350,440,435,1725,3120,650,792.5,1650,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.623765,'amu*angstrom^2'), symmetry=1, barrier=(53.678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.623644,'amu*angstrom^2'), symmetry=1, barrier=(53.6656,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.625226,'amu*angstrom^2'), symmetry=1, barrier=(53.6766,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.61257,0.0445988,-1.22342e-05,-2.00628e-08,1.28136e-11,84409.4,22.2441], Tmin=(100,'K'), Tmax=(907.919,'K')), NASAPolynomial(coeffs=[12.1908,0.0192769,-5.56048e-06,8.55135e-10,-5.54357e-14,81611.4,-32.5986], Tmin=(907.919,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(701.044,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=CC=CCJ) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
)

species(
    label = 'C=C=C(6077)',
    structure = SMILES('C=C=C'),
    E0 = (182.728,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,540,610,2055],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (40.0639,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2746.46,'J/mol'), sigma=(4.78521,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=428.99 K, Pc=56.87 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.36163,0.00777114,2.52438e-05,-3.61894e-08,1.38594e-11,22005.7,7.1245], Tmin=(100,'K'), Tmax=(966.19,'K')), NASAPolynomial(coeffs=[6.46487,0.0106812,-3.73734e-06,6.87011e-10,-4.98581e-14,20670.6,-11.5463], Tmin=(966.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(182.728,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = '[CH]=C=[CH](18734)',
    structure = SMILES('[CH]=C=[CH]'),
    E0 = (491.681,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,239.877,511.233,1743.98,1746.51,1747.6,1753.44],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (38.048,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1737.73,'J/mol'), sigma=(4.1,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.766,0.0170203,-1.57568e-05,7.95984e-09,-1.4265e-12,59188.9,11.2142], Tmin=(100,'K'), Tmax=(1806.04,'K')), NASAPolynomial(coeffs=[4.81405,0.00509933,2.77647e-07,-2.23082e-10,1.96202e-14,59653.5,3.45727], Tmin=(1806.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(491.681,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=C=CJ)"""),
)

species(
    label = '[CH2][C]([CH2])C1[C]=C1(21662)',
    structure = SMILES('[CH2][C]([CH2])C1[C]=C1'),
    E0 = (1008.29,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44701,0.0644249,-0.000101987,9.20436e-08,-3.14104e-11,121353,23.4995], Tmin=(100,'K'), Tmax=(917.641,'K')), NASAPolynomial(coeffs=[4.37614,0.0306842,-1.25514e-05,2.16209e-09,-1.37934e-13,121698,14.4312], Tmin=(917.641,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1008.29,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopropene) + radical(Isobutyl) + radical(cyclopropenyl-vinyl) + radical(Isobutyl) + radical(Tertalkyl)"""),
)

species(
    label = 'H(8)',
    structure = SMILES('[H]'),
    E0 = (211.805,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (1.00794,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1205.6,'J/mol'), sigma=(2.05,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,9.24385e-15,-1.3678e-17,6.66185e-21,-1.00107e-24,25474.2,-0.444973], Tmin=(100,'K'), Tmax=(3459.6,'K')), NASAPolynomial(coeffs=[2.5,9.20456e-12,-3.58608e-15,6.15199e-19,-3.92042e-23,25474.2,-0.444973], Tmin=(3459.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(211.805,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""H""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'C#C[C]=C([CH2])[CH2](21663)',
    structure = SMILES('C#C[C]=C([CH2])[CH2]'),
    E0 = (676.301,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,750,770,3400,2100,350,440,435,1725,2175,525,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100],'cm^-1')),
        HinderedRotor(inertia=(0.112762,'amu*angstrom^2'), symmetry=1, barrier=(79.0845,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.43919,'amu*angstrom^2'), symmetry=1, barrier=(79.0737,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.4403,'amu*angstrom^2'), symmetry=1, barrier=(79.0994,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (77.1039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.62,0.0524854,-5.63586e-05,3.49186e-08,-8.84603e-12,81425.8,17.7705], Tmin=(100,'K'), Tmax=(956.471,'K')), NASAPolynomial(coeffs=[8.94394,0.021856,-8.32304e-06,1.4371e-09,-9.46092e-14,80024.8,-17.2372], Tmin=(956.471,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(676.301,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCtH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(CTCC=CCJ) + radical(Cds_S) + radical(CTCC=CCJ)"""),
)

species(
    label = 'CH2(T)(28)',
    structure = SMILES('[CH2]'),
    E0 = (381.37,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1066.91,2790.99,3622.37],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.01192,-0.000154979,3.26298e-06,-2.40422e-09,5.69497e-13,45867.7,0.5332], Tmin=(100,'K'), Tmax=(1104.58,'K')), NASAPolynomial(coeffs=[3.14983,0.00296674,-9.76056e-07,1.54115e-10,-9.50338e-15,46058.1,4.77808], Tmin=(1104.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH]=C=C[C]=C(19277)',
    structure = SMILES('[CH]=C=C[C]=C'),
    E0 = (587.972,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2950,3100,1380,975,1025,1650,540,610,2055,3010,987.5,1337.5,450,1655,1685,370],'cm^-1')),
        HinderedRotor(inertia=(1.7607,'amu*angstrom^2'), symmetry=1, barrier=(40.4819,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (64.0853,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63528,0.044168,-4.52114e-05,2.40401e-08,-4.84455e-12,70808.6,17.5611], Tmin=(100,'K'), Tmax=(1402.09,'K')), NASAPolynomial(coeffs=[11.4691,0.0097722,-1.62986e-06,9.23713e-11,5.78346e-16,68674.3,-30.9823], Tmin=(1402.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(587.972,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CJC=C) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=C[C]=C([CH2])[CH2](21664)',
    structure = SMILES('[CH]=C[C]=C([CH2])[CH2]'),
    E0 = (701.044,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,1685,370,3010,987.5,1337.5,450,1655,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100],'cm^-1')),
        HinderedRotor(inertia=(0.621136,'amu*angstrom^2'), symmetry=1, barrier=(53.6801,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.62246,'amu*angstrom^2'), symmetry=1, barrier=(53.6952,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.629564,'amu*angstrom^2'), symmetry=1, barrier=(53.6708,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.61259,0.0445986,-1.22336e-05,-2.00637e-08,1.28141e-11,84409.4,22.2441], Tmin=(100,'K'), Tmax=(907.915,'K')), NASAPolynomial(coeffs=[12.1908,0.0192769,-5.56052e-06,8.55145e-10,-5.54365e-14,81611.4,-32.5984], Tmin=(907.915,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(701.044,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(Cds_P) + radical(C=CJC=C) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH]=[C][C]=C([CH2])C(21665)',
    structure = SMILES('[CH]=[C][C]=C([CH2])C'),
    E0 = (781.983,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3120,650,792.5,1650,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(2.23318,'amu*angstrom^2'), symmetry=1, barrier=(51.3451,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.23537,'amu*angstrom^2'), symmetry=1, barrier=(51.3955,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.23516,'amu*angstrom^2'), symmetry=1, barrier=(51.3906,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39604,0.0581521,-5.67597e-05,2.65072e-08,-2.57056e-12,94144.1,21.0741], Tmin=(100,'K'), Tmax=(755.705,'K')), NASAPolynomial(coeffs=[9.82732,0.0236087,-8.21012e-06,1.33516e-09,-8.4593e-14,92581.8,-19.1455], Tmin=(755.705,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(781.983,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=CC=CCJ) + radical(Cds_P) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C([CH2])=[C][C]=C(19938)',
    structure = SMILES('[CH2]C([CH2])=[C][C]=C'),
    E0 = (652.943,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2950,3100,1380,975,1025,1650,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,344.112],'cm^-1')),
        HinderedRotor(inertia=(0.639125,'amu*angstrom^2'), symmetry=1, barrier=(53.7083,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.639213,'amu*angstrom^2'), symmetry=1, barrier=(53.7081,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.639177,'amu*angstrom^2'), symmetry=1, barrier=(53.7083,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14768,0.0530529,-4.37189e-05,2.0291e-08,-3.75323e-12,78641.8,23.1059], Tmin=(100,'K'), Tmax=(1459.24,'K')), NASAPolynomial(coeffs=[11.3261,0.0205991,-5.67812e-06,7.73405e-10,-4.31035e-14,76156.1,-28.1844], Tmin=(1459.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(652.943,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=CC=CCJ) + radical(C=CJC=C) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH]=[C][C]=C([CH2])[CH2](21666)',
    structure = SMILES('[CH]=[C][C]=C([CH2])[CH2]'),
    E0 = (900.039,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,1670,1700,300,440,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100],'cm^-1')),
        HinderedRotor(inertia=(1.71862,'amu*angstrom^2'), symmetry=1, barrier=(53.7158,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.78484,'amu*angstrom^2'), symmetry=1, barrier=(53.675,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.83779,'amu*angstrom^2'), symmetry=1, barrier=(53.6742,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (77.1039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60168,0.0502276,-3.74589e-05,6.6194e-09,4.32853e-12,108339,21.8875], Tmin=(100,'K'), Tmax=(818.743,'K')), NASAPolynomial(coeffs=[10.6982,0.0192977,-5.54693e-06,7.90925e-10,-4.63473e-14,106396,-22.9443], Tmin=(818.743,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(900.039,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=CJC=C) + radical(Cds_P) + radical(C=CC=CCJ) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH]=[C]C1C[C]1[CH2](21667)',
    structure = SMILES('[CH]=[C]C1C[C]1[CH2]'),
    E0 = (956.944,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.79845,0.0429435,-2.9689e-05,1.13377e-08,-1.80342e-12,115178,25.2958], Tmin=(100,'K'), Tmax=(1461.93,'K')), NASAPolynomial(coeffs=[9.66189,0.021428,-7.61301e-06,1.2705e-09,-8.18413e-14,112878,-15.6269], Tmin=(1461.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(956.944,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(Isobutyl) + radical(Cds_S) + radical(Cds_P) + radical(Tertalkyl)"""),
)

species(
    label = '[CH]=C1[CH]C1([CH2])[CH2](21658)',
    structure = SMILES('[CH]=C1[CH]C1([CH2])[CH2]'),
    E0 = (920.691,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.66312,0.0394141,2.57114e-06,-3.42065e-08,1.64357e-11,110828,20.7754], Tmin=(100,'K'), Tmax=(980.909,'K')), NASAPolynomial(coeffs=[14.1709,0.0167535,-6.12013e-06,1.15867e-09,-8.55732e-14,107011,-46.2772], Tmin=(980.909,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(920.691,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsCs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(Neopentyl) + radical(Allyl_S) + radical(Neopentyl) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C1([CH2])[CH][C]=C1(21668)',
    structure = SMILES('[CH2]C1([CH2])[CH][C]=C1'),
    E0 = (907.035,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.05686,0.0295426,2.52825e-05,-5.34587e-08,2.21295e-11,109173,19.864], Tmin=(100,'K'), Tmax=(986.65,'K')), NASAPolynomial(coeffs=[12.5794,0.0189023,-7.22037e-06,1.395e-09,-1.03672e-13,105538,-38.6578], Tmin=(986.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(907.035,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsCs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(cyclobutene-vinyl) + radical(cyclobutene-allyl) + radical(Neopentyl) + radical(Neopentyl)"""),
)

species(
    label = '[CH]=C=C=C([CH2])C(21669)',
    structure = SMILES('[CH]=C=C=C([CH2])C'),
    E0 = (554.202,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26737,0.0571839,-5.52201e-05,2.89993e-08,-6.14086e-12,66756,20.3568], Tmin=(100,'K'), Tmax=(1141.86,'K')), NASAPolynomial(coeffs=[11.5525,0.0211546,-7.89009e-06,1.36592e-09,-9.07513e-14,64407.2,-30.6272], Tmin=(1141.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(554.202,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Cdd-CdsCds) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=C=C[C]1CC1(21670)',
    structure = SMILES('[CH]=C=C[C]1CC1'),
    E0 = (541.359,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28283,0.0223998,4.90804e-05,-8.32727e-08,3.49872e-11,65186.4,17.5341], Tmin=(100,'K'), Tmax=(926.122,'K')), NASAPolynomial(coeffs=[13.0127,0.0154655,-3.51724e-06,5.36516e-10,-3.95877e-14,61508.9,-42.5319], Tmin=(926.122,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(541.359,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + ring(Cyclopropane) + radical(C=C=CJ) + radical(Allyl_T)"""),
)

species(
    label = '[CH]=C1[CH]C(=C)C1(19897)',
    structure = SMILES('[CH]=C1[CH]C(=C)C1'),
    E0 = (544.957,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.58991,0.00539838,0.000111829,-1.55679e-07,6.19498e-11,65617.6,18.851], Tmin=(100,'K'), Tmax=(932.569,'K')), NASAPolynomial(coeffs=[16.877,0.00969419,-5.57109e-07,6.5104e-11,-1.51374e-14,60101.3,-64.3671], Tmin=(932.569,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(544.957,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(C=CCJC=C) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]C=[C][CH2](19274)',
    structure = SMILES('[CH]=[C]C=[C][CH2]'),
    E0 = (859.885,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,1670,1700,300,440,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(2.0544,'amu*angstrom^2'), symmetry=1, barrier=(47.2346,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.05464,'amu*angstrom^2'), symmetry=1, barrier=(47.2402,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (64.0853,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.93557,0.0405981,-3.87335e-05,2.02214e-08,-4.17345e-12,103499,19.2405], Tmin=(100,'K'), Tmax=(1265.13,'K')), NASAPolynomial(coeffs=[9.93801,0.0132358,-3.8481e-06,5.50841e-10,-3.19634e-14,101639,-20.5969], Tmin=(1265.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(859.885,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CC=CCJ) + radical(C=CJC=C) + radical(Cds_P)"""),
)

species(
    label = '[CH]C([CH2])=C(18776)',
    structure = SMILES('[CH]C([CH2])=C'),
    E0 = (489.252,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,350,440,435,1725,3000,3100,440,815,1455,1000,340.08,340.579,341.119],'cm^-1')),
        HinderedRotor(inertia=(0.611693,'amu*angstrom^2'), symmetry=1, barrier=(50.5163,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.611621,'amu*angstrom^2'), symmetry=1, barrier=(50.5272,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (53.0825,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.48527,0.0254697,1.30556e-05,-3.21059e-08,1.33815e-11,58905,14.7489], Tmin=(100,'K'), Tmax=(981.423,'K')), NASAPolynomial(coeffs=[8.34326,0.0203037,-7.64415e-06,1.37964e-09,-9.63641e-14,56854.1,-17.9932], Tmin=(981.423,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(489.252,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(Allyl_P) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[C]=[CH](18830)',
    structure = SMILES('[C]=[CH]'),
    E0 = (847.348,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([562.459,1392.74,3112.83],'cm^-1')),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (25.0293,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.86574,0.0023643,1.31489e-06,-2.64796e-09,1.00932e-12,101918,6.62295], Tmin=(100,'K'), Tmax=(1016.96,'K')), NASAPolynomial(coeffs=[4.17915,0.00251347,-9.43473e-07,1.6872e-10,-1.15831e-14,101782,4.75427], Tmin=(1016.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(847.348,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[CH]C([CH2])=C[C]=[CH](19887)',
    structure = SMILES('[CH]C([CH2])=C[C]=[CH]'),
    E0 = (953.672,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1685,370,350,440,435,1725,3120,650,792.5,1650,3010,987.5,1337.5,450,1655,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.198,'amu*angstrom^2'), symmetry=1, barrier=(50.5364,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.19802,'amu*angstrom^2'), symmetry=1, barrier=(50.5367,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.1981,'amu*angstrom^2'), symmetry=1, barrier=(50.5387,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (77.1039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33019,0.0534685,-3.53423e-05,5.78722e-09,2.75888e-12,114801,22.8075], Tmin=(100,'K'), Tmax=(930.983,'K')), NASAPolynomial(coeffs=[11.325,0.0235591,-8.15178e-06,1.354e-09,-8.87473e-14,112375,-27.7303], Tmin=(930.983,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(953.672,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(C=CJC=C) + radical(Cds_P) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[C]=[C]C=C([CH2])[CH2](21671)',
    structure = SMILES('[C]=[C]C=C([CH2])[CH2]'),
    E0 = (1012.05,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,350,440,435,1725,3010,987.5,1337.5,450,1655,460.913],'cm^-1')),
        HinderedRotor(inertia=(0.357303,'amu*angstrom^2'), symmetry=1, barrier=(55.0896,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.363325,'amu*angstrom^2'), symmetry=1, barrier=(55.1357,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.36112,'amu*angstrom^2'), symmetry=1, barrier=(55.0833,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (77.1039,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36111,0.0508615,-4.47237e-05,2.15625e-08,-4.15269e-12,121823,23.1206], Tmin=(100,'K'), Tmax=(1329.46,'K')), NASAPolynomial(coeffs=[11.6217,0.0180679,-5.55461e-06,8.33303e-10,-5.01259e-14,119264,-28.6636], Tmin=(1329.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1012.05,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(C=CC=CCJ) + radical(C=CJC=C) + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[CH][C]=[CH](21256)',
    structure = SMILES('[CH][C]=[CH]'),
    E0 = (861.746,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.1891,'amu*angstrom^2'), symmetry=1, barrier=(50.3317,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (38.048,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.18317,0.0164338,-7.13252e-06,1.19383e-09,-3.27944e-14,103675,12.0918], Tmin=(100,'K'), Tmax=(1799.19,'K')), NASAPolynomial(coeffs=[6.32962,0.0112581,-4.33439e-06,7.19107e-10,-4.49321e-14,102248,-5.75439], Tmin=(1799.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(861.746,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2][C]=C(6078)',
    structure = SMILES('[CH2][C]=C'),
    E0 = (395.465,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,1685,370],'cm^-1')),
        HinderedRotor(inertia=(0.243012,'amu*angstrom^2'), symmetry=1, barrier=(30.4931,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (40.0639,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.28772,0.0102496,1.79964e-05,-2.86352e-08,1.11956e-11,47593.9,10.4766], Tmin=(100,'K'), Tmax=(969.996,'K')), NASAPolynomial(coeffs=[6.37267,0.0109726,-3.91218e-06,7.11399e-10,-5.07602e-14,46362.9,-7.57277], Tmin=(969.996,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(395.465,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH][C](C)C=C=[CH](19893)',
    structure = SMILES('[CH][C](C)C=C=[CH]'),
    E0 = (826.977,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,540,610,2055,3010,987.5,1337.5,450,1655,360,370,350,363.861,363.869,3828.89],'cm^-1')),
        HinderedRotor(inertia=(0.0422098,'amu*angstrom^2'), symmetry=1, barrier=(84.0671,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.894785,'amu*angstrom^2'), symmetry=1, barrier=(84.0672,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.894776,'amu*angstrom^2'), symmetry=1, barrier=(84.0672,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59501,0.0467612,-2.52228e-05,-2.84185e-09,5.5575e-12,99554.6,21.4529], Tmin=(100,'K'), Tmax=(937.475,'K')), NASAPolynomial(coeffs=[11.3962,0.0200355,-6.61135e-06,1.09739e-09,-7.29451e-14,97053.6,-28.7369], Tmin=(937.475,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(826.977,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CCJ2_triplet) + radical(C=C=CJ) + radical(Allyl_T)"""),
)

species(
    label = '[CH]C([CH2])=CC=[CH](19885)',
    structure = SMILES('[CH]C([CH2])=CC=[CH]'),
    E0 = (754.677,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3120,650,792.5,1650,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,196.348,196.487,196.558],'cm^-1')),
        HinderedRotor(inertia=(1.83328,'amu*angstrom^2'), symmetry=1, barrier=(50.2305,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.83418,'amu*angstrom^2'), symmetry=1, barrier=(50.2304,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.83143,'amu*angstrom^2'), symmetry=1, barrier=(50.2326,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39638,0.0470718,-6.79632e-06,-2.63294e-08,1.41795e-11,90869.6,22.9727], Tmin=(100,'K'), Tmax=(955.331,'K')), NASAPolynomial(coeffs=[12.8884,0.0234311,-8.10912e-06,1.40591e-09,-9.68773e-14,87556.9,-37.7903], Tmin=(955.331,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(754.677,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(AllylJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[CH][C]=C[C]1CC1(21672)',
    structure = SMILES('[CH][C]=C[C]1CC1'),
    E0 = (818.804,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.20321,0.0262219,4.18197e-05,-6.89443e-08,2.72817e-11,98556,19.7219], Tmin=(100,'K'), Tmax=(966.213,'K')), NASAPolynomial(coeffs=[10.1114,0.0263868,-9.51824e-06,1.7233e-09,-1.22468e-13,95491.9,-26.107], Tmin=(966.213,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(818.804,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopropane) + radical(AllylJ2_triplet) + radical(Allyl_T) + radical(Cds_S)"""),
)

species(
    label = '[CH]C1=C[C]([CH2])C1(19935)',
    structure = SMILES('[CH]C1=C[C]([CH2])C1'),
    E0 = (781.408,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.19791,0.0249446,5.08976e-05,-8.29823e-08,3.38632e-11,94060,20.362], Tmin=(100,'K'), Tmax=(934.201,'K')), NASAPolynomial(coeffs=[11.1082,0.0240103,-7.36034e-06,1.23652e-09,-8.65644e-14,90771.2,-30.7108], Tmin=(934.201,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(781.408,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(Allyl_T) + radical(Isobutyl) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH][C]=[C]C([CH2])[CH2](21443)',
    structure = SMILES('[CH][C]=[C]C([CH2])[CH2]'),
    E0 = (1172.23,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,633.003,633.003,633.005,633.005],'cm^-1')),
        HinderedRotor(inertia=(0.198658,'amu*angstrom^2'), symmetry=1, barrier=(56.4868,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.198658,'amu*angstrom^2'), symmetry=1, barrier=(56.4868,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.198659,'amu*angstrom^2'), symmetry=1, barrier=(56.4868,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.19866,'amu*angstrom^2'), symmetry=1, barrier=(56.4868,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32288,0.060526,-6.68314e-05,4.68359e-08,-1.36415e-11,141082,28.4193], Tmin=(100,'K'), Tmax=(923.89,'K')), NASAPolynomial(coeffs=[7.16022,0.0309136,-1.1708e-05,1.97569e-09,-1.26833e-13,140188,1.72195], Tmin=(923.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1172.23,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cds_S) + radical(Isobutyl) + radical(AllylJ2_triplet) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = '[CH][C]=CC[C]=C(19650)',
    structure = SMILES('[CH][C]=CC[C]=C'),
    E0 = (899.092,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1670,1700,300,440,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,373.303,380.079,382.075,383.044,383.267],'cm^-1')),
        HinderedRotor(inertia=(0.487812,'amu*angstrom^2'), symmetry=1, barrier=(51.6967,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.488,'amu*angstrom^2'), symmetry=1, barrier=(51.691,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.498368,'amu*angstrom^2'), symmetry=1, barrier=(51.7361,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3289.69,'J/mol'), sigma=(5.72141,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=513.84 K, Pc=39.86 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.88044,0.0486308,-3.16094e-05,1.10825e-08,-1.71752e-12,108210,24.2874], Tmin=(100,'K'), Tmax=(1376.51,'K')), NASAPolynomial(coeffs=[7.40561,0.0325751,-1.41134e-05,2.60888e-09,-1.78548e-13,106689,-4.13394], Tmin=(1376.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(899.092,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(AllylJ2_triplet) + radical(Cds_S)"""),
)

species(
    label = '[CH]C(=C)C[C]=[CH](18821)',
    structure = SMILES('[CH]C(=C)C[C]=[CH]'),
    E0 = (905.317,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3120,650,792.5,1650,2950,3100,1380,975,1025,1650,351.212,351.261,351.665,351.728],'cm^-1')),
        HinderedRotor(inertia=(0.581927,'amu*angstrom^2'), symmetry=1, barrier=(51.055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.582794,'amu*angstrom^2'), symmetry=1, barrier=(51.0423,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.581489,'amu*angstrom^2'), symmetry=1, barrier=(51.0473,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3226.68,'J/mol'), sigma=(5.67268,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=504.00 K, Pc=40.11 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46847,0.0539599,-4.0095e-05,1.60924e-08,-2.72865e-12,108977,24.6512], Tmin=(100,'K'), Tmax=(1348.57,'K')), NASAPolynomial(coeffs=[10.0073,0.0286331,-1.19248e-05,2.16668e-09,-1.47115e-13,106674,-19.0976], Tmin=(1348.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(905.317,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(AllylJ2_triplet) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]C([CH2])=C[C]=C(18814)',
    structure = SMILES('[CH]C([CH2])=C[C]=C'),
    E0 = (706.576,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,350,440,435,1725,1685,370,3010,987.5,1337.5,450,1655,213.788,218.272,220.282,221.947],'cm^-1')),
        HinderedRotor(inertia=(1.47498,'amu*angstrom^2'), symmetry=1, barrier=(50.7819,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.49216,'amu*angstrom^2'), symmetry=1, barrier=(50.8696,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.53172,'amu*angstrom^2'), symmetry=1, barrier=(50.8545,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3312.87,'J/mol'), sigma=(5.74688,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=517.46 K, Pc=39.61 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39389,0.0500704,-1.91482e-05,-1.09538e-08,8.38261e-12,85081.9,22.1753], Tmin=(100,'K'), Tmax=(944.304,'K')), NASAPolynomial(coeffs=[11.1102,0.0263476,-9.15979e-06,1.54654e-09,-1.03094e-13,82469.5,-28.2594], Tmin=(944.304,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(706.576,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(AllylJ2_triplet) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2][C]1C=[C][CH]C1(21673)',
    structure = SMILES('[CH2][C]1C=[C][CH]C1'),
    E0 = (712.12,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.75288,0.0112327,7.3533e-05,-1.01804e-07,3.96495e-11,85708,17.2359], Tmin=(100,'K'), Tmax=(940.166,'K')), NASAPolynomial(coeffs=[10.4304,0.0197448,-5.74374e-06,9.95147e-10,-7.32702e-14,82444.5,-29.0085], Tmin=(940.166,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(712.12,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopentene) + radical(Isobutyl) + radical(cyclopentene-vinyl) + radical(Allyl_T) + radical(cyclopentene-allyl)"""),
)

species(
    label = '[CH]=[C]C1CC1=C(21656)',
    structure = SMILES('[CH]=[C]C1CC1=C'),
    E0 = (743.584,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47833,0.0489614,-3.76714e-05,1.50812e-08,-2.43949e-12,89528.7,19.0576], Tmin=(100,'K'), Tmax=(1464.97,'K')), NASAPolynomial(coeffs=[12.3201,0.0193587,-7.3609e-06,1.28778e-09,-8.56202e-14,86352.1,-37.3875], Tmin=(1464.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(743.584,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(Cds_P) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C=[C]C([CH2])[CH2](21440)',
    structure = SMILES('[CH]=C=[C]C([CH2])[CH2]'),
    E0 = (894.788,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,1685,370,540,610,2055,1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100],'cm^-1')),
        HinderedRotor(inertia=(3.05512,'amu*angstrom^2'), symmetry=1, barrier=(70.2433,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00230006,'amu*angstrom^2'), symmetry=1, barrier=(9.68046,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.05756,'amu*angstrom^2'), symmetry=1, barrier=(70.2992,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23552,0.0587671,-6.72418e-05,4.31429e-08,-1.0793e-11,107719,26.8243], Tmin=(100,'K'), Tmax=(1105.19,'K')), NASAPolynomial(coeffs=[10.3789,0.0194281,-5.3719e-06,7.08084e-10,-3.7148e-14,106080,-16.476], Tmin=(1105.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(894.788,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cdd-Cd)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Isobutyl) + radical(C=C=CJ) + radical(Cds_S) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C#CC([CH2])=C(20384)',
    structure = SMILES('[CH2]C#CC([CH2])=C'),
    E0 = (495.977,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,350,440,435,1725,2100,2250,500,550,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,402.479],'cm^-1')),
        HinderedRotor(inertia=(0.0355994,'amu*angstrom^2'), symmetry=1, barrier=(25.836,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0356068,'amu*angstrom^2'), symmetry=1, barrier=(25.8373,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.66726,'amu*angstrom^2'), symmetry=1, barrier=(84.3175,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.7757,0.0405469,-1.02263e-05,-1.50206e-08,8.45145e-12,59739.6,21.2046], Tmin=(100,'K'), Tmax=(1007.86,'K')), NASAPolynomial(coeffs=[11.3952,0.0206698,-7.88038e-06,1.44424e-09,-1.01686e-13,56871.1,-29.8906], Tmin=(1007.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(495.977,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-CtHHH) + group(Cds-CdsCtCs) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(Allyl_P) + radical(Propargyl)"""),
)

species(
    label = '[CH2][C][CH2](10272)',
    structure = SMILES('[CH2][C][CH2]'),
    E0 = (738.184,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,349.736],'cm^-1')),
        HinderedRotor(inertia=(0.00138263,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0013768,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (40.0639,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.07912,0.0169806,-2.30739e-06,-7.22553e-09,3.53372e-12,88819,13.4505], Tmin=(100,'K'), Tmax=(1024.31,'K')), NASAPolynomial(coeffs=[6.32415,0.0113321,-4.3212e-06,7.79353e-10,-5.38459e-14,87785.8,-4.0815], Tmin=(1024.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(738.184,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJ) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[C]=[C]CC([CH2])=C(19668)',
    structure = SMILES('[C]=[C]CC([CH2])=C'),
    E0 = (997.138,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,350,440,435,1725,3000,3100,440,815,1455,1000,392.967,798.686],'cm^-1')),
        HinderedRotor(inertia=(0.024156,'amu*angstrom^2'), symmetry=1, barrier=(10.9229,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.477392,'amu*angstrom^2'), symmetry=1, barrier=(10.9762,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.2249,'amu*angstrom^2'), symmetry=1, barrier=(28.1629,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60508,0.0531797,-4.81003e-05,2.36592e-08,-4.83642e-12,120014,22.6182], Tmin=(100,'K'), Tmax=(1152.18,'K')), NASAPolynomial(coeffs=[9.87995,0.0244519,-1.06999e-05,2.0187e-09,-1.40848e-13,118107,-18.4754], Tmin=(1152.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(997.138,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CdCdJ2_triplet) + radical(Cds_S) + radical(Allyl_P)"""),
)

species(
    label = '[C]#C[CH]C([CH2])[CH2](21674)',
    structure = SMILES('[C]#C[CH]C([CH2])[CH2]'),
    E0 = (986.498,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,1380,1390,370,380,2900,435,2175,525,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.0041152,'amu*angstrom^2'), symmetry=1, barrier=(5.0308,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00411531,'amu*angstrom^2'), symmetry=1, barrier=(5.0265,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.74004,'amu*angstrom^2'), symmetry=1, barrier=(62.9989,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0515551,'amu*angstrom^2'), symmetry=1, barrier=(62.9992,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73665,0.053266,-3.95242e-05,-1.48769e-08,2.82404e-11,118727,23.7812], Tmin=(100,'K'), Tmax=(615.708,'K')), NASAPolynomial(coeffs=[10.4228,0.0175306,-2.88174e-06,3.74275e-11,1.93283e-14,117265,-17.0971], Tmin=(615.708,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(986.498,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Isobutyl) + radical(Isobutyl) + radical(Acetyl) + radical(Sec_Propargyl)"""),
)

species(
    label = '[C]=[C]C=C([CH2])C(21675)',
    structure = SMILES('[C]=[C]C=C([CH2])C'),
    E0 = (893.994,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,1685,370,3010,987.5,1337.5,450,1655,180.001],'cm^-1')),
        HinderedRotor(inertia=(1.07302,'amu*angstrom^2'), symmetry=1, barrier=(24.6708,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.28031,'amu*angstrom^2'), symmetry=1, barrier=(75.4208,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.84751,'amu*angstrom^2'), symmetry=1, barrier=(75.427,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40493,0.0559529,-5.47087e-05,3.04206e-08,-6.93847e-12,107617,21.4058], Tmin=(100,'K'), Tmax=(1055.48,'K')), NASAPolynomial(coeffs=[9.88747,0.0238063,-9.02348e-06,1.56476e-09,-1.03703e-13,105826,-19.9756], Tmin=(1055.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(893.994,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(CdCdJ2_triplet) + radical(C=CJC=C)"""),
)

species(
    label = '[C]#C(5143)',
    structure = SMILES('[C]#C'),
    E0 = (551.936,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (25.0293,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.03853,0.0115449,-2.13263e-05,1.81932e-08,-5.41587e-12,66398,5.96675], Tmin=(100,'K'), Tmax=(1076.58,'K')), NASAPolynomial(coeffs=[4.00845,0.00206817,6.04935e-08,-1.17707e-10,1.2928e-14,66529.5,2.79656], Tmin=(1076.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(551.936,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(62.3585,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Ct-CtH) + group(Ct-CtH) + radical(Acetyl)"""),
)

species(
    label = '[CH][C]([CH2])[CH2](21676)',
    structure = SMILES('[CH][C]([CH2])[CH2]'),
    E0 = (890.457,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([360,370,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1461.46,1461.87,1462.89],'cm^-1')),
        HinderedRotor(inertia=(0.0040705,'amu*angstrom^2'), symmetry=1, barrier=(6.18251,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.269078,'amu*angstrom^2'), symmetry=1, barrier=(6.18664,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00408396,'amu*angstrom^2'), symmetry=1, barrier=(6.18869,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (53.0825,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.47236,0.0388202,-5.79863e-05,5.22017e-08,-1.76428e-11,107147,19.8163], Tmin=(100,'K'), Tmax=(949.729,'K')), NASAPolynomial(coeffs=[3.20586,0.0214581,-8.02212e-06,1.30527e-09,-7.98424e-14,107652,19.7043], Tmin=(949.729,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(890.457,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(CCJ2_triplet) + radical(Isobutyl) + radical(Tertalkyl)"""),
)

species(
    label = 'N2',
    structure = SMILES('N#N'),
    E0 = (-8.64289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0135,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(810.913,'J/mol'), sigma=(3.621,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.76,'angstroms^3'), rotrelaxcollnum=4.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53101,-0.000123661,-5.02999e-07,2.43531e-09,-1.40881e-12,-1046.98,2.96747], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.95258,0.0013969,-4.92632e-07,7.8601e-11,-4.60755e-15,-923.949,5.87189], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-8.64289,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""N2""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'Ne',
    structure = SMILES('[Ne]'),
    E0 = (-6.19738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (20.1797,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1235.53,'J/mol'), sigma=(3.758e-10,'m'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with fixed Lennard Jones Parameters. This is the fallback method! Try improving transport databases!"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,3.35532], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,3.35532], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-6.19738,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""Ne""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'He',
    structure = SMILES('[He]'),
    E0 = (-6.19738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (4.0026,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(84.8076,'J/mol'), sigma=(2.576,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,0.928724], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,0.928724], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-6.19738,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""He""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'Ar',
    structure = SMILES('[Ar]'),
    E0 = (-6.19738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (39.348,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1134.93,'J/mol'), sigma=(3.33,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,4.37967], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,4.37967], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-6.19738,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""Ar""", comment="""Thermo library: primaryThermoLibrary"""),
)

transitionState(
    label = 'TS1',
    E0 = (701.044,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (1010.2,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (895.88,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (983.995,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (863.823,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (937.94,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (923.369,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (1111.84,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (956.944,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (920.691,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (907.035,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (764.444,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (709.133,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (709.328,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (1275.57,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (1370.92,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (1165.48,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (1223.85,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (1076.5,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (909.831,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (956.292,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (1158.11,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (1257.21,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (818.804,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (839.099,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (1195.09,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (1069.13,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (1097.45,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (897.785,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (733.397,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (743.584,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (1072.47,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (764.444,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (1264.18,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (1181.31,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (1138.97,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (935.565,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (1476.71,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=[C]C=C([CH2])[CH2](19662)'],
    products = ['C=C=C(6077)', '[CH]=C=[CH](18734)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=[C]C=C([CH2])[CH2](19662)'],
    products = ['[CH2][C]([CH2])C1[C]=C1(21662)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.54e+10,'s^-1'), n=0.69, Ea=(309.157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(8)', 'C#C[C]=C([CH2])[CH2](21663)'],
    products = ['[CH]=[C]C=C([CH2])[CH2](19662)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CH2(T)(28)', '[CH]=C=C[C]=C(19277)'],
    products = ['[CH]=[C]C=C([CH2])[CH2](19662)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.70446,'m^3/(mol*s)'), n=2.07639, Ea=(14.6531,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ca_Cds-HH;YJ] for rate rule [Ca_Cds-HH;CH2_triplet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]=C[C]=C([CH2])[CH2](21664)'],
    products = ['[CH]=[C]C=C([CH2])[CH2](19662)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]=[C][C]=C([CH2])C(21665)'],
    products = ['[CH]=[C]C=C([CH2])[CH2](19662)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(5.17353e+06,'s^-1'), n=1.89718, Ea=(155.956,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=[C]C=C([CH2])[CH2](19662)'],
    products = ['[CH2]C([CH2])=[C][C]=C(19938)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(191.5,'s^-1'), n=3.05, Ea=(222.325,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleH;Cd_H_out_doubleC] for rate rule [R3HJ;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(8)', '[CH]=[C][C]=C([CH2])[CH2](21666)'],
    products = ['[CH]=[C]C=C([CH2])[CH2](19662)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=[C]C=C([CH2])[CH2](19662)'],
    products = ['[CH]=[C]C1C[C]1[CH2](21667)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4.07074e+13,'s^-1'), n=-0.296394, Ea=(255.9,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 135 used for R3_D;doublebond_intra;radadd_intra_cs2H
Exact match found for rate rule [R3_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 254.8 to 255.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=[C]C=C([CH2])[CH2](19662)'],
    products = ['[CH]=C1[CH]C1([CH2])[CH2](21658)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.42595e+10,'s^-1'), n=0.7335, Ea=(219.648,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3_D;doublebond_intra_pri;radadd_intra]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 218.7 to 219.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=[C]C=C([CH2])[CH2](19662)'],
    products = ['[CH2]C1([CH2])[CH][C]=C1(21668)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(205.991,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_linear;doublebond_intra_pri;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 203.7 to 206.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=[C]C=C([CH2])[CH2](19662)'],
    products = ['[CH]=C=C=C([CH2])C(21669)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=[C]C=C([CH2])[CH2](19662)'],
    products = ['[CH]=C=C[C]1CC1(21670)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(7.38971e+10,'s^-1'), n=0.0476667, Ea=(8.08907,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_2H;Cpri_rad_out_2H] for rate rule [R3_SS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=[C]C=C([CH2])[CH2](19662)'],
    products = ['[CH]=C1[CH]C(=C)C1(19897)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SDS;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['CH2(T)(28)', '[CH]=[C]C=[C][CH2](19274)'],
    products = ['[CH]=[C]C=C([CH2])[CH2](19662)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]C([CH2])=C(18776)', '[C]=[CH](18830)'],
    products = ['[CH]=[C]C=C([CH2])[CH2](19662)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(8)', '[CH]C([CH2])=C[C]=[CH](19887)'],
    products = ['[CH]=[C]C=C([CH2])[CH2](19662)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(8)', '[C]=[C]C=C([CH2])[CH2](21671)'],
    products = ['[CH]=[C]C=C([CH2])[CH2](19662)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C=C=C(6077)', '[CH][C]=[CH](21256)'],
    products = ['[CH]=[C]C=C([CH2])[CH2](19662)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(0.00173894,'m^3/(mol*s)'), n=2.67356, Ea=(32.0272,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ca_Cds-HH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2][C]=C(6078)', '[CH]=C=[CH](18734)'],
    products = ['[CH]=[C]C=C([CH2])[CH2](19662)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.04713,'m^3/(mol*s)'), n=2.10494, Ea=(22.6844,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH][C](C)C=C=[CH](19893)'],
    products = ['[CH]=[C]C=C([CH2])[CH2](19662)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(720,'s^-1'), n=2.932, Ea=(129.315,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 57 used for R3H_DS;Cd_rad_out_singleH;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]C([CH2])=CC=[CH](19885)'],
    products = ['[CH]=[C]C=C([CH2])[CH2](19662)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.13764e+12,'s^-1'), n=1.09983, Ea=(403.433,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSD;Cd_rad_out_single;XH_out] for rate rule [R4H_DSD;Cd_rad_out_singleH;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2][C]=C(6078)', '[CH][C]=[CH](21256)'],
    products = ['[CH]=[C]C=C([CH2])[CH2](19662)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(3.9578e+07,'m^3/(mol*s)'), n=-0.126319, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Y_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination
Ea raised from -15.6 to -15.6 kJ/mol.
Ea raised from -15.6 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH]=[C]C=C([CH2])[CH2](19662)'],
    products = ['[CH][C]=C[C]1CC1(21672)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.03537e+13,'s^-1'), n=-0.296394, Ea=(117.76,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra_cs2H] for rate rule [R3_D;doublebond_intra_secDe_2H;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 115.5 to 117.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH]=[C]C=C([CH2])[CH2](19662)'],
    products = ['[CH]C1=C[C]([CH2])C1(19935)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(3.906e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra;radadd_intra]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH][C]=[C]C([CH2])[CH2](21443)'],
    products = ['[CH]=[C]C=C([CH2])[CH2](19662)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH][C]=CC[C]=C(19650)'],
    products = ['[CH]=[C]C=C([CH2])[CH2](19662)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]C(=C)C[C]=[CH](18821)'],
    products = ['[CH]=[C]C=C([CH2])[CH2](19662)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(13437.7,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH]C([CH2])=C[C]=C(18814)'],
    products = ['[CH]=[C]C=C([CH2])[CH2](19662)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.456e+11,'s^-1'), n=0.86, Ea=(191.209,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cd_H_out_singleH] for rate rule [R5Hall;Cd_rad_out_singleH;Cd_H_out_singleH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH]=[C]C=C([CH2])[CH2](19662)'],
    products = ['[CH2][C]1C=[C][CH]C1(21673)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.82795e+11,'s^-1'), n=0.0396934, Ea=(32.3534,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_linear;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH]=[C]C=C([CH2])[CH2](19662)'],
    products = ['[CH]=[C]C1CC1=C(21656)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.18842e+14,'s^-1'), n=0.0123667, Ea=(42.5404,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_out;Cpri_rad_out_2H] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination
Ea raised from 42.5 to 42.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH]=C=[C]C([CH2])[CH2](21440)'],
    products = ['[CH]=[C]C=C([CH2])[CH2](19662)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(3.9024e+09,'s^-1'), n=1.09833, Ea=(177.686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_double;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH]=[C]C=C([CH2])[CH2](19662)'],
    products = ['[CH2]C#CC([CH2])=C(20384)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2][C][CH2](10272)', '[CH]=C=[CH](18734)'],
    products = ['[CH]=[C]C=C([CH2])[CH2](19662)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(2.29708e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[C]=[C]CC([CH2])=C(19668)'],
    products = ['[CH]=[C]C=C([CH2])[CH2](19662)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(2.30972e+07,'s^-1'), n=1.70216, Ea=(184.173,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_MS;Y_rad_out;XH_out] for rate rule [R3H_TS;Ct_rad_out;XH_out]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[C]#C[CH]C([CH2])[CH2](21674)'],
    products = ['[CH]=[C]C=C([CH2])[CH2](19662)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(322699,'s^-1'), n=1.96429, Ea=(152.474,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Y_rad_out;XH_out] for rate rule [R4HJ_2;Ct_rad_out;XH_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[C]=[C]C=C([CH2])C(21675)'],
    products = ['[CH]=[C]C=C([CH2])[CH2](19662)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(505536,'s^-1'), n=1.7378, Ea=(41.5716,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Y_rad_out;Cs_H_out_2H] for rate rule [R5Hall;Ct_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[C]#C(5143)', '[CH][C]([CH2])[CH2](21676)'],
    products = ['[CH]=[C]C=C([CH2])[CH2](19662)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Ct_rad;Birad]
Euclidian distance = 1.0
family: Birad_R_Recombination"""),
)

network(
    label = '4533',
    isomers = [
        '[CH]=[C]C=C([CH2])[CH2](19662)',
    ],
    reactants = [
        ('C=C=C(6077)', '[CH]=C=[CH](18734)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '4533',
    Tmin = (1200,'K'),
    Tmax = (1500,'K'),
    Tcount = 10,
    Tlist = ([1201.48,1213.22,1236.21,1269.31,1310.55,1356.92,1404.16,1447.02,1479.84,1497.7],'K'),
    Pmin = (1,'atm'),
    Pmax = (10,'atm'),
    Pcount = 10,
    Plist = ([1.02771,1.14872,1.41959,1.89986,2.67608,3.83649,5.40396,7.23219,8.93758,9.98989],'bar'),
    maximumGrainSize = (0.5,'kcal/mol'),
    minimumGrainCount = 250,
    method = 'modified strong collision',
    interpolationModel = ('Chebyshev', 6, 4),
    activeKRotor = True,
    activeJRotor = True,
    rmgmode = True,
)

