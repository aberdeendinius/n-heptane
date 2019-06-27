species(
    label = 'C[CH]OC=C=C[O](22453)',
    structure = SMILES('C[CH]OC=C=C[O]'),
    E0 = (81.6691,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,540,610,2055,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.16343,'amu*angstrom^2'), symmetry=1, barrier=(26.7497,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.16054,'amu*angstrom^2'), symmetry=1, barrier=(26.6831,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.16003,'amu*angstrom^2'), symmetry=1, barrier=(26.6714,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0573684,0.0710985,-3.48287e-05,-3.07917e-08,2.41938e-11,9985.56,24.851], Tmin=(100,'K'), Tmax=(913.113,'K')), NASAPolynomial(coeffs=[27.4917,-0.00103565,3.91744e-06,-8.54088e-10,5.54982e-14,2930.58,-116.636], Tmin=(913.113,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(81.6691,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(CCsJOC(O)) + radical(C=COJ)"""),
)

species(
    label = 'CC=O(606)',
    structure = SMILES('CC=O'),
    E0 = (-177.906,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,1427.17,1427.17,1427.17,1427.17,3755.47],'cm^-1')),
        HinderedRotor(inertia=(0.717734,'amu*angstrom^2'), symmetry=1, barrier=(16.5021,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.70079,0.000387835,3.86929e-05,-4.52447e-08,1.58859e-11,-21380.9,9.13562], Tmin=(100,'K'), Tmax=(984.198,'K')), NASAPolynomial(coeffs=[4.58889,0.0128894,-4.91502e-06,9.26508e-10,-6.71011e-14,-22336,0.901072], Tmin=(984.198,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-177.906,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsH)"""),
)

species(
    label = 'C#CC=O(21959)',
    structure = SMILES('C#CC=O'),
    E0 = (84.2941,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,750,770,3400,2100,346.083,346.148,1254.41,1569.49,1569.49],'cm^-1')),
        HinderedRotor(inertia=(0.289847,'amu*angstrom^2'), symmetry=1, barrier=(24.6472,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (54.0474,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3153.34,'J/mol'), sigma=(5.09602,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=492.54 K, Pc=54.07 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.00307,0.020875,-1.61807e-05,6.30357e-09,-1.00066e-12,10174.9,9.76778], Tmin=(100,'K'), Tmax=(1462.04,'K')), NASAPolynomial(coeffs=[7.33456,0.00902438,-4.02236e-06,7.59526e-10,-5.26541e-14,8908.32,-12.7744], Tmin=(1462.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(84.2941,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-O2d(Cds-Cds)H) + group(Ct-CtCs) + group(Ct-CtH)"""),
)

species(
    label = 'CC1OC1[C]=C[O](23242)',
    structure = SMILES('CC1OC1[C]=C[O]'),
    E0 = (134.324,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16251,0.0454109,1.84374e-05,-7.51422e-08,3.90321e-11,16273.8,23.5259], Tmin=(100,'K'), Tmax=(877.638,'K')), NASAPolynomial(coeffs=[20.3392,0.00619878,3.0943e-06,-9.24126e-10,6.92528e-14,11051.9,-77.0607], Tmin=(877.638,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(134.324,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Ethylene_oxide) + radical(C=COJ) + radical(Cds_S)"""),
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
    label = 'C=COC=C=C[O](22353)',
    structure = SMILES('C=COC=C=C[O]'),
    E0 = (43.0625,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,540,610,2055,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.36209,'amu*angstrom^2'), symmetry=1, barrier=(31.317,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.36961,'amu*angstrom^2'), symmetry=1, barrier=(31.4901,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.77056,0.0537876,-7.49554e-06,-4.6763e-08,2.69792e-11,5311.51,24.7606], Tmin=(100,'K'), Tmax=(921.336,'K')), NASAPolynomial(coeffs=[22.7797,0.00320453,1.64259e-06,-3.98035e-10,2.31728e-14,-652.726,-89.9759], Tmin=(921.336,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(43.0625,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=COJ)"""),
)

species(
    label = 'C[CH]OC=C=C=O(23243)',
    structure = SMILES('C[CH]OC=C=C=O'),
    E0 = (125.766,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,540,610,2055,3010,987.5,1337.5,450,1655,2120,512.5,787.5,3025,407.5,1350,352.5,180,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.320845,'amu*angstrom^2'), symmetry=1, barrier=(7.37686,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.324108,'amu*angstrom^2'), symmetry=1, barrier=(7.45188,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.321999,'amu*angstrom^2'), symmetry=1, barrier=(7.40338,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22132,0.0475551,-1.66371e-05,-2.84822e-08,1.89051e-11,15238.8,7.4626], Tmin=(100,'K'), Tmax=(922.514,'K')), NASAPolynomial(coeffs=[20.2509,0.000490464,2.25283e-06,-4.81186e-10,2.91693e-14,10219.4,-90.9845], Tmin=(922.514,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(125.766,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(CCsJOC(O))"""),
)

species(
    label = '[CH]=C=C[O](8556)',
    structure = SMILES('[CH]=C=C[O]'),
    E0 = (269.875,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,3010,987.5,1337.5,450,1655,540,610,2055],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0474,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.7401,0.0176211,1.62543e-05,-4.58526e-08,2.291e-11,32513.4,12.64], Tmin=(100,'K'), Tmax=(883.628,'K')), NASAPolynomial(coeffs=[13.5244,-0.00323064,4.17673e-06,-9.22668e-10,6.45049e-14,29515.7,-44.2318], Tmin=(883.628,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(269.875,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C=CJ) + radical(C=COJ)"""),
)

species(
    label = '[CH2]COC=C=C[O](23244)',
    structure = SMILES('[CH2]COC=C=C[O]'),
    E0 = (99.3308,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,540,610,2055,3000,3100,440,815,1455,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.19274,'amu*angstrom^2'), symmetry=1, barrier=(27.4234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.19617,'amu*angstrom^2'), symmetry=1, barrier=(27.5023,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.19589,'amu*angstrom^2'), symmetry=1, barrier=(27.4959,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.156282,0.0669155,-2.78351e-05,-3.39989e-08,2.42921e-11,12101.6,25.2438], Tmin=(100,'K'), Tmax=(916.623,'K')), NASAPolynomial(coeffs=[25.7985,0.00167101,2.58716e-06,-5.97797e-10,3.77459e-14,5440.78,-106.923], Tmin=(916.623,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(99.3308,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(CJCO) + radical(C=COJ)"""),
)

species(
    label = 'C[CH]OC=C=[C]O(23245)',
    structure = SMILES('C[CH]OC=C=[C]O'),
    E0 = (179.951,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,540,610,2055,3025,407.5,1350,352.5,1685,370,3010,987.5,1337.5,450,1655,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.00534,'amu*angstrom^2'), symmetry=1, barrier=(23.1148,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00351,'amu*angstrom^2'), symmetry=1, barrier=(23.0726,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00636,'amu*angstrom^2'), symmetry=1, barrier=(23.1383,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00402,'amu*angstrom^2'), symmetry=1, barrier=(23.0844,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.27975,0.0918499,-0.000107204,5.66199e-08,-1.08726e-11,21854,30.9571], Tmin=(100,'K'), Tmax=(1506.15,'K')), NASAPolynomial(coeffs=[26.5948,-0.00106614,4.14341e-06,-9.91641e-10,7.21821e-14,15599.7,-107.826], Tmin=(1506.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(179.951,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=CJO) + radical(CCsJOC(O))"""),
)

species(
    label = 'CCO[C]=C=C[O](23246)',
    structure = SMILES('CCO[C]=C=C[O]'),
    E0 = (127.486,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2850,1437.5,1250,1305,750,350,540,610,2055,3010,987.5,1337.5,450,1655,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.03345,'amu*angstrom^2'), symmetry=1, barrier=(23.7611,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03401,'amu*angstrom^2'), symmetry=1, barrier=(23.774,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03558,'amu*angstrom^2'), symmetry=1, barrier=(23.81,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.516069,0.0635236,-3.51727e-05,-1.39826e-08,1.39893e-11,15470.5,26.4466], Tmin=(100,'K'), Tmax=(933.537,'K')), NASAPolynomial(coeffs=[20.8557,0.00955772,-1.78254e-06,2.51122e-10,-1.98896e-14,10226.9,-78.0266], Tmin=(933.537,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(127.486,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = 'C[CH]O[C]=C=CO(23247)',
    structure = SMILES('C[CH]O[C]=C=CO'),
    E0 = (179.951,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,540,610,2055,3025,407.5,1350,352.5,1685,370,3010,987.5,1337.5,450,1655,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.00551,'amu*angstrom^2'), symmetry=1, barrier=(23.1186,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00532,'amu*angstrom^2'), symmetry=1, barrier=(23.1142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00217,'amu*angstrom^2'), symmetry=1, barrier=(23.0418,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00345,'amu*angstrom^2'), symmetry=1, barrier=(23.0712,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.27978,0.0918502,-0.000107204,5.66208e-08,-1.08728e-11,21854,30.9572], Tmin=(100,'K'), Tmax=(1506.12,'K')), NASAPolynomial(coeffs=[26.5956,-0.00106737,4.14404e-06,-9.91777e-10,7.21926e-14,15599.2,-107.831], Tmin=(1506.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(179.951,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(CCsJOC(O)) + radical(C=CJO)"""),
)

species(
    label = 'CCO[CH][C]=C=O(23248)',
    structure = SMILES('CCO[CH][C]=C=O'),
    E0 = (94.4529,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.898954,0.0783476,-0.000127875,1.19582e-07,-4.32383e-11,11462,25.6538], Tmin=(100,'K'), Tmax=(851.793,'K')), NASAPolynomial(coeffs=[4.74596,0.0382227,-1.83691e-05,3.47195e-09,-2.36717e-13,11606.8,12.4088], Tmin=(851.793,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(94.4529,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsOsHH) + group(Cs-(Cds-Cdd-O2d)OsHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + radical(Cds_S) + radical(C=CCJ(O)C)"""),
)

species(
    label = 'C=CO[CH][C]=CO(23249)',
    structure = SMILES('C=CO[CH][C]=CO'),
    E0 = (64.0403,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.370538,0.0622611,-1.92257e-05,-4.05937e-08,2.62296e-11,7849.34,26.0826], Tmin=(100,'K'), Tmax=(914.392,'K')), NASAPolynomial(coeffs=[24.7091,0.00243488,2.40126e-06,-5.76481e-10,3.6777e-14,1448.42,-99.8211], Tmin=(914.392,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(64.0403,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CCJ(O)C) + radical(Cds_S)"""),
)

species(
    label = 'C[CH][O](605)',
    structure = SMILES('C[CH][O]'),
    E0 = (149.432,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,2066.51],'cm^-1')),
        HinderedRotor(inertia=(0.362113,'amu*angstrom^2'), symmetry=1, barrier=(8.32568,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.20363,0.021847,-3.14755e-05,3.43227e-08,-1.42322e-11,17997,11.0861], Tmin=(100,'K'), Tmax=(846.374,'K')), NASAPolynomial(coeffs=[1.2024,0.020386,-9.53523e-06,1.79858e-09,-1.23081e-13,18726.8,22.7175], Tmin=(846.374,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(149.432,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCsJOH) + radical(CCOJ)"""),
)

species(
    label = 'C=CO[CH][C]=C[O](22369)',
    structure = SMILES('C=CO[CH][C]=C[O]'),
    E0 = (205.503,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,243.301,243.301,243.301,243.302,243.302],'cm^-1')),
        HinderedRotor(inertia=(0.710148,'amu*angstrom^2'), symmetry=1, barrier=(29.8308,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.710148,'amu*angstrom^2'), symmetry=1, barrier=(29.8308,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.710148,'amu*angstrom^2'), symmetry=1, barrier=(29.8308,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.784875,0.055271,-1.58329e-05,-3.4209e-08,2.13262e-11,24846.3,26.0047], Tmin=(100,'K'), Tmax=(932.038,'K')), NASAPolynomial(coeffs=[21.381,0.00624748,-2.9389e-07,-4.90047e-12,-4.20798e-15,19297.1,-81.0833], Tmin=(932.038,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(205.503,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CCJ(O)C) + radical(Cds_S) + radical(C=COJ)"""),
)

species(
    label = 'C[CH]O[C]=C=C[O](23250)',
    structure = SMILES('C[CH]O[C]=C=C[O]'),
    E0 = (321.413,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,540,610,2055,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.05516,'amu*angstrom^2'), symmetry=1, barrier=(24.2602,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05503,'amu*angstrom^2'), symmetry=1, barrier=(24.2572,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05628,'amu*angstrom^2'), symmetry=1, barrier=(24.2858,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.613875,0.0818389,-9.29393e-05,4.83655e-08,-9.2553e-12,38840.1,29.98], Tmin=(100,'K'), Tmax=(1481,'K')), NASAPolynomial(coeffs=[24.0794,0.00162522,1.99707e-06,-5.33952e-10,3.97126e-14,33008.6,-93.8427], Tmin=(1481,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(321.413,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=CJO) + radical(CCsJOC(O)) + radical(C=COJ)"""),
)

species(
    label = 'C[CH]O[CH][C]=C=O(23251)',
    structure = SMILES('C[CH]O[CH][C]=C=O'),
    E0 = (274.909,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2120,512.5,787.5,3000,3050,390,425,1340,1360,335,370,272.122,275.388,275.751,3298.75],'cm^-1')),
        HinderedRotor(inertia=(0.177095,'amu*angstrom^2'), symmetry=1, barrier=(9.27467,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.909711,'amu*angstrom^2'), symmetry=1, barrier=(48.9868,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.171083,'amu*angstrom^2'), symmetry=1, barrier=(9.27843,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.922246,'amu*angstrom^2'), symmetry=1, barrier=(48.964,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.574483,0.0856194,-0.000147791,1.33929e-07,-4.61765e-11,33177.4,26.2649], Tmin=(100,'K'), Tmax=(874.611,'K')), NASAPolynomial(coeffs=[7.62179,0.0306745,-1.46018e-05,2.71199e-09,-1.81592e-13,32813.5,-1.82355], Tmin=(874.611,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(274.909,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsOsHH) + group(Cs-(Cds-Cdd-O2d)OsHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + radical(CCsJOCs) + radical(Cds_S) + radical(C=CCJ(O)C)"""),
)

species(
    label = 'C[CH]OC=C1[CH]O1(23252)',
    structure = SMILES('C[CH]OC=C1[CH]O1'),
    E0 = (107.851,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.320908,0.0504,4.07463e-05,-1.21258e-07,5.93764e-11,13132.5,22.0749], Tmin=(100,'K'), Tmax=(911.906,'K')), NASAPolynomial(coeffs=[32.8116,-0.011561,1.01579e-05,-2.02306e-09,1.30289e-13,3857.38,-150.043], Tmin=(911.906,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(107.851,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(methyleneoxirane) + radical(C=CCJO) + radical(CCsJOC(O))"""),
)

species(
    label = 'CC1O[CH]C1=C[O](23253)',
    structure = SMILES('CC1O[CH]C1=C[O]'),
    E0 = (10.853,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80186,0.0266915,6.22894e-05,-1.0456e-07,4.27641e-11,1404.09,23.4947], Tmin=(100,'K'), Tmax=(959.581,'K')), NASAPolynomial(coeffs=[17.5437,0.0153313,-4.76986e-06,9.55846e-10,-7.82798e-14,-4115.12,-64.8177], Tmin=(959.581,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(10.853,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + ring(Cyclobutane) + radical(C=CCJ(O)C) + radical(C=COJ)"""),
)

species(
    label = 'C[CH]OC1[C]=CO1(23254)',
    structure = SMILES('C[CH]OC1[C]=CO1'),
    E0 = (199.393,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.547071,0.0603315,-2.27515e-05,-2.81874e-08,1.90135e-11,24120,22.9401], Tmin=(100,'K'), Tmax=(944.183,'K')), NASAPolynomial(coeffs=[21.9995,0.00819305,-1.47216e-06,2.4802e-10,-2.29162e-14,18342,-88.4688], Tmin=(944.183,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(199.393,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(Cds_S) + radical(CCsJOCs)"""),
)

species(
    label = 'CC1OC=[C]C1[O](23255)',
    structure = SMILES('CC1OC=[C]C1[O]'),
    E0 = (130.449,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56277,0.0306296,5.99214e-05,-1.16137e-07,5.21409e-11,15798.7,19.1507], Tmin=(100,'K'), Tmax=(906.923,'K')), NASAPolynomial(coeffs=[21.2398,0.0043962,3.1592e-06,-7.91964e-10,5.16026e-14,9739.28,-87.5865], Tmin=(906.923,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(130.449,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(2,3-Dihydrofuran) + radical(CC(C)OJ) + radical(Cds_S)"""),
)

species(
    label = 'CCOC=C=C=O(23256)',
    structure = SMILES('CCOC=C=C=O'),
    E0 = (-68.1617,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.58517,0.0382077,1.00197e-05,-5.0738e-08,2.51131e-11,-8097.31,6.68257], Tmin=(100,'K'), Tmax=(935.976,'K')), NASAPolynomial(coeffs=[17.769,0.00693534,-5.87535e-07,6.93654e-11,-1.0254e-14,-12786.6,-79.1906], Tmin=(935.976,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-68.1617,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C=COC=C=CO(23257)',
    structure = SMILES('C=COC=C=CO'),
    E0 = (-98.4001,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.354641,0.0607922,-1.09142e-05,-5.31617e-08,3.19157e-11,-11685.4,24.8443], Tmin=(100,'K'), Tmax=(907.324,'K')), NASAPolynomial(coeffs=[26.1326,-0.000649898,4.3617e-06,-9.7525e-10,6.46234e-14,-18511.9,-108.854], Tmin=(907.324,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-98.4001,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds)"""),
)

species(
    label = 'C[CH]O[C]=C[CH][O](23258)',
    structure = SMILES('C[CH]O[C]=C[CH][O]'),
    E0 = (402.003,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3000,3050,390,425,1340,1360,335,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.440849,0.0688378,-6.69043e-05,3.22197e-08,-6.03982e-12,48486,31.0862], Tmin=(100,'K'), Tmax=(1306.93,'K')), NASAPolynomial(coeffs=[17.8657,0.015507,-5.69503e-06,9.96717e-10,-6.72395e-14,43931.3,-57.6431], Tmin=(1306.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(402.003,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=CCJO) + radical(CCsJOC(O)) + radical(C=CJO) + radical(CCOJ)"""),
)

species(
    label = 'C[CH]O[CH]C=[C][O](23259)',
    structure = SMILES('C[CH]O[CH]C=[C][O]'),
    E0 = (299.995,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3000,3050,390,425,1340,1360,335,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1271,0.0665334,-7.40337e-05,4.65575e-08,-1.21351e-11,36181.7,29.1911], Tmin=(100,'K'), Tmax=(920.672,'K')), NASAPolynomial(coeffs=[9.72019,0.0292,-1.32096e-05,2.51492e-09,-1.75951e-13,34599.4,-11.5555], Tmin=(920.672,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(299.995,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CCsJOCs) + radical(C=CCJ(O)C) + radical(C=CJO)"""),
)

species(
    label = '[CH]C(32)',
    structure = SMILES('[CH]C'),
    E0 = (351.472,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,431.535,1804.51],'cm^-1')),
        HinderedRotor(inertia=(0.000906356,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.73285,-0.000244773,3.59198e-05,-4.44289e-08,1.65887e-11,42287.5,7.078], Tmin=(100,'K'), Tmax=(940.479,'K')), NASAPolynomial(coeffs=[5.42969,0.0081677,-2.4253e-06,4.22642e-10,-3.09417e-14,41277.1,-4.6789], Tmin=(940.479,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(351.472,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ2_triplet)"""),
)

species(
    label = '[O]C=C=C[O](22349)',
    structure = SMILES('[O]C=C=C[O]'),
    E0 = (48.0679,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,540,610,2055,180,180],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0468,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.10387,0.0254095,2.29542e-05,-6.61705e-08,3.24893e-11,5864.8,15.5654], Tmin=(100,'K'), Tmax=(907.068,'K')), NASAPolynomial(coeffs=[19.4725,-0.00784929,6.29342e-06,-1.25746e-09,8.23938e-14,931.198,-76.3609], Tmin=(907.068,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(48.0679,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = '[CH3](11)',
    structure = SMILES('[CH3]'),
    E0 = (135.382,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([570.572,1408.13,1408.49,4000,4000,4000],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (15.0345,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.91547,0.00184154,3.48742e-06,-3.32748e-09,8.49957e-13,16285.6,0.351741], Tmin=(100,'K'), Tmax=(1337.63,'K')), NASAPolynomial(coeffs=[3.54146,0.00476787,-1.82148e-06,3.28877e-10,-2.22546e-14,16224,1.66035], Tmin=(1337.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(135.382,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo library: primaryThermoLibrary + radical(CH3)"""),
)

species(
    label = '[CH]OC=C=C[O](23145)',
    structure = SMILES('[CH]OC=C=C[O]'),
    E0 = (386.278,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,540,610,2055,180,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.0316,'amu*angstrom^2'), symmetry=1, barrier=(23.7185,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03202,'amu*angstrom^2'), symmetry=1, barrier=(23.7281,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.947942,0.0502343,-1.58912e-05,-3.98424e-08,2.5868e-11,46584.2,19.7696], Tmin=(100,'K'), Tmax=(908.33,'K')), NASAPolynomial(coeffs=[24.7588,-0.00908885,6.88222e-06,-1.37011e-09,9.0156e-14,40380.2,-103.154], Tmin=(908.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(386.278,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-OsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(CH2_triplet)"""),
)

species(
    label = 'O(T)(63)',
    structure = SMILES('[O]'),
    E0 = (243.034,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (15.9994,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(665.16,'J/mol'), sigma=(2.75,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,9.24385e-15,-1.3678e-17,6.66185e-21,-1.00107e-24,29230.2,4.09104], Tmin=(100,'K'), Tmax=(3459.6,'K')), NASAPolynomial(coeffs=[2.5,9.20456e-12,-3.58608e-15,6.15199e-19,-3.92042e-23,29230.2,4.09104], Tmin=(3459.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.034,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""O(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH]=C=CO[CH]C(19117)',
    structure = SMILES('[CH]=C=CO[CH]C'),
    E0 = (303.476,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3120,650,792.5,1650,2750,2800,2850,1350,1500,750,1050,1375,1000,540,610,2055,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.12591,'amu*angstrom^2'), symmetry=1, barrier=(25.8868,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12676,'amu*angstrom^2'), symmetry=1, barrier=(25.9065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12744,'amu*angstrom^2'), symmetry=1, barrier=(25.922,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.399911,0.0747907,-8.15405e-05,4.15088e-08,-7.72286e-12,36677,25.4427], Tmin=(100,'K'), Tmax=(1564.34,'K')), NASAPolynomial(coeffs=[21.208,0.00361349,1.98028e-06,-5.92768e-10,4.55074e-14,31865.3,-82.2433], Tmin=(1564.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(303.476,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CCsJOC(O)) + radical(C=C=CJ)"""),
)

species(
    label = 'C[C]OC=C=C[O](23260)',
    structure = SMILES('C[C]OC=C=C[O]'),
    E0 = (348.73,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,540,610,2055,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.20567,'amu*angstrom^2'), symmetry=1, barrier=(27.7206,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.19831,'amu*angstrom^2'), symmetry=1, barrier=(27.5514,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.20748,'amu*angstrom^2'), symmetry=1, barrier=(27.7623,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.108105,0.0675208,-3.26621e-05,-3.07466e-08,2.36232e-11,42099.5,23.5152], Tmin=(100,'K'), Tmax=(918.037,'K')), NASAPolynomial(coeffs=[27.217,-0.00283231,4.24711e-06,-8.76209e-10,5.55687e-14,35109.4,-115.914], Tmin=(918.037,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(348.73,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(CH2_triplet)"""),
)

species(
    label = 'C[CH]OC#CC=O(23261)',
    structure = SMILES('C[CH]OC#CC=O'),
    E0 = (74.3147,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2782.5,750,1395,475,1775,1000,2100,2250,500,550,3025,407.5,1350,352.5,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.0576,'amu*angstrom^2'), symmetry=1, barrier=(24.3163,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05492,'amu*angstrom^2'), symmetry=1, barrier=(24.2547,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05769,'amu*angstrom^2'), symmetry=1, barrier=(24.3184,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05413,'amu*angstrom^2'), symmetry=1, barrier=(24.2365,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.460424,0.0692428,-6.65439e-05,2.61786e-08,-2.35032e-12,9073.06,23.1177], Tmin=(100,'K'), Tmax=(999.316,'K')), NASAPolynomial(coeffs=[18.9735,0.0113735,-4.04786e-06,7.41998e-10,-5.35866e-14,4562.41,-70.2401], Tmin=(999.316,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(74.3147,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-Cds)H) + group(Ct-CtCs) + group(Ct-CtOs) + radical(CCsJOCs)"""),
)

species(
    label = 'C[CH]O[C]=CC=O(23262)',
    structure = SMILES('C[CH]O[C]=CC=O'),
    E0 = (124.677,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,2750,2800,2850,1350,1500,750,1050,1375,1000,2782.5,750,1395,475,1775,1000,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.699303,'amu*angstrom^2'), symmetry=1, barrier=(16.0783,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.698235,'amu*angstrom^2'), symmetry=1, barrier=(16.0538,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.69876,'amu*angstrom^2'), symmetry=1, barrier=(16.0659,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.698111,'amu*angstrom^2'), symmetry=1, barrier=(16.0509,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.650895,0.0727541,-7.78432e-05,4.14491e-08,-8.73937e-12,15116.5,26.3881], Tmin=(100,'K'), Tmax=(1149.05,'K')), NASAPolynomial(coeffs=[15.6477,0.0205489,-9.69433e-06,1.91044e-09,-1.37057e-13,11670.1,-48.0468], Tmin=(1149.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(124.677,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cd-Cd(CO)H) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + radical(C=CJO) + radical(CCsJOC(O))"""),
)

species(
    label = 'C[CH]O[CH]C=C=O(23263)',
    structure = SMILES('C[CH]O[CH]C=C=O'),
    E0 = (37.0671,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.631457,0.0807233,-0.00012191,1.03183e-07,-3.45689e-11,4573.08,25.5877], Tmin=(100,'K'), Tmax=(844.245,'K')), NASAPolynomial(coeffs=[8.69498,0.0313784,-1.44445e-05,2.69187e-09,-1.82689e-13,3608.57,-9.59762], Tmin=(844.245,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(37.0671,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsOsHH) + group(Cs-(Cds-Cdd-O2d)OsHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + radical(C=CCJ(O)C) + radical(CCsJOCs)"""),
)

species(
    label = 'C=CO[CH]C=C[O](14357)',
    structure = SMILES('C=CO[CH]C=C[O]'),
    E0 = (-32.3389,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,218.463,220.112,223.343,223.704,223.709],'cm^-1')),
        HinderedRotor(inertia=(0.835872,'amu*angstrom^2'), symmetry=1, barrier=(30.4554,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.866639,'amu*angstrom^2'), symmetry=1, barrier=(30.4333,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.84839,'amu*angstrom^2'), symmetry=1, barrier=(30.473,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.812861,0.0507271,8.783e-06,-6.33059e-08,3.22481e-11,-3756.85,25.4306], Tmin=(100,'K'), Tmax=(931.876,'K')), NASAPolynomial(coeffs=[22.5921,0.00671119,4.43052e-09,-5.87494e-11,-2.48098e-15,-9963.91,-89.6298], Tmin=(931.876,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-32.3389,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=CCJ(O)C)"""),
)

species(
    label = 'CC1O[CH][C]=CO1(23264)',
    structure = SMILES('CC1O[CH][C]=CO1'),
    E0 = (26.9426,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73403,0.0279206,6.31493e-05,-1.13526e-07,4.93187e-11,3342.42,17.9356], Tmin=(100,'K'), Tmax=(917.312,'K')), NASAPolynomial(coeffs=[18.9955,0.00952964,2.13286e-07,-1.91669e-10,8.79529e-15,-2217.49,-76.8953], Tmin=(917.312,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(26.9426,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(24dihydro13dioxin) + radical(C=CCJ(O)C) + radical(Cds_S)"""),
)

species(
    label = 'C=COC=CC=O(14375)',
    structure = SMILES('C=COC=CC=O'),
    E0 = (-153.674,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.915281,0.0578288,-3.79768e-05,5.08222e-09,2.52457e-12,-18363.1,25.1928], Tmin=(100,'K'), Tmax=(1081.27,'K')), NASAPolynomial(coeffs=[16.1716,0.0188351,-8.08304e-06,1.57196e-09,-1.13688e-13,-22682.2,-54.318], Tmin=(1081.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-153.674,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cd-Cd(CO)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH)"""),
)

species(
    label = 'CC1OC=C1C=O(23265)',
    structure = SMILES('CC1OC=C1C=O'),
    E0 = (-168.303,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23425,0.037336,4.39984e-05,-9.48455e-08,4.12565e-11,-20121.3,20.3292], Tmin=(100,'K'), Tmax=(960.513,'K')), NASAPolynomial(coeffs=[22.2042,0.00844486,-2.14231e-06,5.20012e-10,-5.10944e-14,-26845.3,-94.0262], Tmin=(960.513,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-168.303,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cd-CdCs(CO)) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + ring(Cyclobutene)"""),
)

species(
    label = '[CH]=O(373)',
    structure = SMILES('[CH]=O'),
    E0 = (33.3613,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([880.033,2045.55,4000],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (29.018,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.62,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.35602,-0.00347087,1.25664e-05,-9.99488e-09,2.27888e-12,3995.77,2.75113], Tmin=(100,'K'), Tmax=(1565.72,'K')), NASAPolynomial(coeffs=[4.61867,0.00504456,-4.39241e-06,9.73282e-10,-7.07436e-14,2787.5,-2.22961], Tmin=(1565.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(33.3613,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-OdHH) + radical(HCdsJO)"""),
)

species(
    label = '[C]=CO[CH]C(5275)',
    structure = SMILES('[C]=CO[CH]C'),
    E0 = (566.522,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.914459,'amu*angstrom^2'), symmetry=1, barrier=(21.0252,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.91447,'amu*angstrom^2'), symmetry=1, barrier=(21.0255,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.914253,'amu*angstrom^2'), symmetry=1, barrier=(21.0205,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (69.0819,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.490641,0.0610275,-6.7446e-05,3.43148e-08,-6.41432e-12,68277.1,21.1379], Tmin=(100,'K'), Tmax=(1526.39,'K')), NASAPolynomial(coeffs=[18.9588,0.00181762,1.36622e-06,-3.80652e-10,2.83888e-14,63898.9,-71.6444], Tmin=(1526.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(566.522,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(CCsJOC(O)) + radical(CdCdJ2_triplet)"""),
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
    E0 = (81.6691,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (172.167,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (261.327,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (340.492,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (127.855,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (257.695,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (357.637,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (268.905,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (495.113,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (108.551,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (188.57,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (419.307,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (417.308,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (533.218,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (486.714,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (268.868,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (207.026,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (207.327,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (139.827,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (90.0371,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (106.642,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (424.304,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (322.857,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (404.706,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (555.975,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (710.357,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (560.535,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (301.229,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (321.136,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (244.448,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (114.709,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (108.028,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (120.894,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (89.5768,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (634.199,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C[CH]OC=C=C[O](22453)'],
    products = ['CC=O(606)', 'C#CC=O(21959)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C[CH]OC=C=C[O](22453)'],
    products = ['CC1OC1[C]=C[O](23242)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(8.15385e+14,'s^-1'), n=-0.537569, Ea=(90.4981,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_csHNd]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(8)', 'C=COC=C=C[O](22353)'],
    products = ['C[CH]OC=C=C[O](22453)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.67e+12,'cm^3/(mol*s)'), n=0.1, Ea=(6.4601,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 2816 used for Cds-HH_Cds-OsH;HJ
Exact match found for rate rule [Cds-HH_Cds-OsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(8)', 'C[CH]OC=C=C=O(23243)'],
    products = ['C[CH]OC=C=C[O](22453)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(7048.5,'m^3/(mol*s)'), n=1.43107, Ea=(2.92155,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CC=O(606)', '[CH]=C=C[O](8556)'],
    products = ['C[CH]OC=C=C[O](22453)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4e+09,'cm^3/(mol*s)'), n=1.39, Ea=(35.8862,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Od_CO-CsH;YJ] for rate rule [Od_CO-CsH;CdsJ=Cdd]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]COC=C=C[O](23244)'],
    products = ['C[CH]OC=C=C[O](22453)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.7e+13,'s^-1','+|-',2), n=-0.1, Ea=(158.364,'kJ/mol'), T0=(1,'K'), Tmin=(700,'K'), Tmax=(1800,'K'), comment="""From training reaction 347 used for R2H_S;C_rad_out_2H;Cs_H_out_H/NonDeO
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/NonDeO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C[CH]OC=C=[C]O(23245)'],
    products = ['C[CH]OC=C=C[O](22453)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(3.9024e+09,'s^-1'), n=1.09833, Ea=(177.686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_double;XH_out] for rate rule [R2H_S;Cd_rad_out_double;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CCO[C]=C=C[O](23246)'],
    products = ['C[CH]OC=C=C[O](22453)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(3.32e+09,'s^-1'), n=0.99, Ea=(141.419,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Cd_rad_out_double;Cs_H_out_H/NonDeC] for rate rule [R3H_SS_O;Cd_rad_out_double;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C[CH]O[C]=C=CO(23247)'],
    products = ['C[CH]OC=C=C[O](22453)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3293.16,'s^-1'), n=2.5965, Ea=(315.162,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_MMS;Cd_rad_out_single;XH_out] for rate rule [R4H_MMS;Cd_rad_out_singleNd;O_H_out]
Euclidian distance = 1.41421356237
family: intra_H_migration"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C[CH]OC=C=C[O](22453)'],
    products = ['CCO[CH][C]=C=O(23248)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(0.0378492,'s^-1'), n=3.26, Ea=(26.8822,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5H;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C[CH]OC=C=C[O](22453)'],
    products = ['C=CO[CH][C]=CO(23249)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(5.85e+08,'s^-1'), n=0, Ea=(106.901,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R7Hall;O_rad_out;Cs_H_out_2H] for rate rule [R7HJ_5;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C[CH][O](605)', '[CH]=C=C[O](8556)'],
    products = ['C[CH]OC=C=C[O](22453)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(40577.9,'m^3/(mol*s)'), n=0.702818, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_allenic;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -11.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(8)', 'C=CO[CH][C]=C[O](22369)'],
    products = ['C[CH]OC=C=C[O](22453)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(8)', 'C[CH]O[C]=C=C[O](23250)'],
    products = ['C[CH]OC=C=C[O](22453)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.048e+13,'cm^3/(mol*s)'), n=0.206, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 44 used for Cd_allenic;H_rad
Exact match found for rate rule [Cd_allenic;H_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -0.7 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(8)', 'C[CH]O[CH][C]=C=O(23251)'],
    products = ['C[CH]OC=C=C[O](22453)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(4.34078e+06,'m^3/(mol*s)'), n=0.278577, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Y_rad]
Euclidian distance = 0
family: R_Recombination
Ea raised from -1.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction16',
    reactants = ['C[CH]OC=C=C[O](22453)'],
    products = ['C[CH]OC=C1[CH]O1(23252)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(6.59543e+11,'s^-1'), n=0.400725, Ea=(187.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra] for rate rule [R3_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C[CH]OC=C=C[O](22453)'],
    products = ['CC1O[CH]C1=C[O](23253)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(5.16207e+08,'s^-1'), n=0.911389, Ea=(125.357,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_csHCs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C[CH]OC=C=C[O](22453)'],
    products = ['C[CH]OC1[C]=CO1(23254)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.47276e+11,'s^-1'), n=0.183155, Ea=(125.658,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_CdCdd;radadd_intra] + [R4;doublebond_intra;radadd_intra] for rate rule [R4;doublebond_intra_CdCdd;radadd_intra_O]
Euclidian distance = 1.41421356237
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C[CH]OC=C=C[O](22453)'],
    products = ['CC1OC=[C]C1[O](23255)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4.64e+06,'s^-1'), n=1.15, Ea=(58.1576,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;doublebond_intra;radadd_intra_csHCs] for rate rule [R5_linear;doublebond_intra_CdCdd;radadd_intra_csHCs]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C[CH]OC=C=C[O](22453)'],
    products = ['CCOC=C=C=O(23256)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(3.21e+09,'s^-1'), n=0.137, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5;Y_rad_NDe;XH_Rrad] for rate rule [R5radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C[CH]OC=C=C[O](22453)'],
    products = ['C=COC=C=CO(23257)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;Y_rad;XH_Rrad] for rate rule [R7radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C[CH]O[C]=C[CH][O](23258)'],
    products = ['C[CH]OC=C=C[O](22453)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C[CH]O[CH]C=[C][O](23259)'],
    products = ['C[CH]OC=C=C[O](22453)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH]C(32)', '[O]C=C=C[O](22349)'],
    products = ['C[CH]OC=C=C[O](22453)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(87.1677,'m^3/(mol*s)'), n=1.88017, Ea=(5.1666,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [O_sec_rad;Birad] for rate rule [O_rad/OneDe;Birad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH3](11)', '[CH]OC=C=C[O](23145)'],
    products = ['C[CH]OC=C=C[O](22453)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_methyl;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['O(T)(63)', '[CH]=C=CO[CH]C(19117)'],
    products = ['C[CH]OC=C=C[O](22453)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(93609.6,'m^3/(mol*s)'), n=1.13083, Ea=(163.847,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['H(8)', 'C[C]OC=C=C[O](23260)'],
    products = ['C[CH]OC=C=C[O](22453)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['H(8)', 'C[CH]OC#CC=O(23261)'],
    products = ['C[CH]OC=C=C[O](22453)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct_Ct;HJ] for rate rule [Ct-O_Ct;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C[CH]O[C]=CC=O(23262)'],
    products = ['C[CH]OC=C=C[O](22453)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.75e+11,'s^-1'), n=0.633, Ea=(196.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_singleNd;Cd_H_out_singleDe]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C[CH]OC=C=C[O](22453)'],
    products = ['C[CH]O[CH]C=C=O(23263)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C[CH]OC=C=C[O](22453)'],
    products = ['C=CO[CH]C=C[O](14357)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_single;Cs_H_out_2H] for rate rule [R5HJ_3;Cd_rad_out_singleDe;Cs_H_out_2H]
Euclidian distance = 1.41421356237
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C[CH]OC=C=C[O](22453)'],
    products = ['CC1O[CH][C]=CO1(23264)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(487000,'s^-1'), n=1.17, Ea=(26.3592,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;multiplebond_intra;radadd_intra_csHCs] for rate rule [R6_linear;carbonyl_intra_H;radadd_intra_csHCs]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C[CH]OC=C=C[O](22453)'],
    products = ['C=COC=CC=O(14375)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(5.55988e+09,'s^-1'), n=0.137, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_De;XH_Rrad] for rate rule [R5radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C[CH]OC=C=C[O](22453)'],
    products = ['CC1OC=C1C=O(23265)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_1H;Ypri_rad_out] + [R4;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_H/NonDeC;CdsinglepriDe_rad_out]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH]=O(373)', '[C]=CO[CH]C(5275)'],
    products = ['C[CH]OC=C=C[O](22453)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.14854e+06,'m^3/(mol*s)'), n=0.575199, Ea=(34.3157,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

network(
    label = '4749',
    isomers = [
        'C[CH]OC=C=C[O](22453)',
    ],
    reactants = [
        ('CC=O(606)', 'C#CC=O(21959)'),
    ],
    bathGas = {
        'N2': 0.25,
        'Ne': 0.25,
        'He': 0.25,
        'Ar': 0.25,
    },
)

pressureDependence(
    label = '4749',
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

