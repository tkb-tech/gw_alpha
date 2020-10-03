import gwdet

print('please input m1')
m1 = input()
print('please input m2')
m2 = input()
print('please input z')
z = input()

p=gwdet.detectability()
#m1=10. # Component mass in Msun
#m2=10. # Component mass in Msun
#z=0.1  # Redshift

print(p(m1,m2,z))  # Fraction of detectabile sources
