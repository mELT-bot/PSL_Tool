from dataclasses import dataclass




@dataclass
class Material:
    name:str
    young_modulus_mpa:int
    poisson_ratio:float


#https://www.sonelastic.com/en/fundamentals/tables-of-materials-properties/non-ferrous-metals.html
Aluminum=Material('Aluminum Alloy 1100',69000,0.33)
Titanium=Material('Titanium Alloy Ti- 6A1-4V',114000,0.34)
Gold=Material('Gold (Pure)',171000,0.39)
Led=Material('Chemical Led',13500,0.44)

#https://www.engineersedge.com/manufacturing_spec/average_properties_structural_materials.htm
CastIron=Material('Cast Iron',99973,.21)

#https://www.sonelastic.com/en/fundamentals/tables-of-materials-properties/polymers.html
Nylon=Material('Nylon 6.6',3000,0.39)
Polycarbonate=Material('Polycarbonate',2380,0.36)

#https://www.sonelastic.com/en/fundamentals/tables-of-materials-properties/ferrous-metals.html
StainlessSteel=Material('Stainless Alloy 304',193000,0.30)

#https://www.sonelastic.com/en/fundamentals/tables-of-materials-properties/woods.html
Wood=Material('Eucalipto Grandis (Eucalyptus grandis)',12800,0.30)

#Deniz
Concrete=Material('Concrete',29000,0.2)