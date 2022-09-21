

            

            

m1 = Mirror(0, 3, 0.005)
m2 = Mirror(0.009, 30, 0.01)
m3 = Mirror(0.267, 0, 0.25)


Res = Resonator(3, m1, m2, m3)

print(Res.is_consistent())
print(Res.is_g_stable_sagittal())
print(Res.is_g_stable_tangential())
print(Res.get_length())
print(Res.elems[-1].angle)

Res.system_scheme()

Res.realign()


print(Res.is_consistent())
print(Res.is_g_stable_sagittal())
print(Res.is_g_stable_tangential())
print(Res.get_length())
print(Res.waist_search(1064 * 10 **-9))

Res.system_scheme()
