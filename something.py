#import bpy  # Импорт основных модулей Blender Python API
from math import *
def f():
    print("Вы в первом методе. Результаты:")
    def Z(x, y):
        return round(pow((x * cos(30) + (y - 1) * sin(30)), 2) / 4 + pow(((y - 1) * cos(30) + x * sin(30)), 2) / 36, 3)
        # return round(100*pow((y-pow(x,2)),2)+pow((1-x),2), 3)

    def C(x1, x2):
        return 0.5 * (x2 + x1)

    def New(x, c):
        return 2 * c - x

    def R(x, c):
        return -c + 2 * x

    def S(x, c):
        return 0.5 * (x + c)

    def half(x, y):
        return x + 0.5 * (x - y)

    V = ((-3.3, 3.5, Z(-3.3, 3.5)), (-3.3, 3.4, Z(-3.3, 3.4)), (-3.2, 3.5, Z(-3.2, 3.5)))
    Xh = ()
    Xl = ()
    Xg = ()
    Xn = ()
    Xc = ()
    Xr = ()
    Xs = ()
    Q = 1
    N = 0
    VERTEXES = ()
    while Q >= 0.001:
        #for i in range(0, 3):
            #bpy.ops.mesh.primitive_plane_add(radius=0.005, location=(V[i]))
        if V[0][2] > V[1][2] and V[0][2] > V[2][2]:
            Xh = V[0]
            if V[1][2] > V[2][2]:
                Xg = V[1]
                Xl = V[2]
            else:
                Xg = V[2]
                Xl = V[1]
        elif V[1][2] > V[0][2] and V[1][2] > V[2][2]:
            Xh = V[1]
            if V[0][2] > V[2][2]:
                Xg = V[0]
                Xl = V[2]
            else:
                Xg = V[2]
                Xl = V[0]
        else:
            Xh = V[2]
            if V[0][2] > V[1][2]:
                Xg = V[0]
                Xl = V[1]
            else:
                Xg = V[1]
                Xl = V[0]
        VERTEXES += Xh,
        VERTEXES += Xg,
        VERTEXES += Xl,
        Xc = (C(Xg[0], Xl[0]), C(Xg[1], Xl[1]), Z(C(Xg[0], Xl[0]), C(Xg[1], Xl[1])))
        X0 = (New(Xh[0], Xc[0]), New(Xh[1], Xc[1]), Z(New(Xh[0], Xc[0]), New(Xh[1], Xc[1])))
        # вывод вершин мб в полигона в будущем
        #bpy.ops.mesh.primitive_cylinder_add(vertices=6, radius=0.005, depth=0.1, location=(Xc))
        #bpy.ops.mesh.primitive_plane_add(radius=0.005, location=(X0))

        if X0[2] < Xl[2]:
            Xr = R(X0[0], Xc[0]), R(X0[1], Xc[1]), Z(R(X0[0], Xc[0]), R(X0[1], Xc[1]))
            if Xr[2] < Xl[2]:
                Xh = Xr
            else:
                Xh = X0
            #bpy.ops.mesh.primitive_cylinder_add(vertices=3, radius=0.01, depth=0.1, location=(Xr))
        elif Xl[2] < X0[2] <= Xg[2]:
            Xh = X0
        elif X0[2] > Xg[2]:
            if X0[2] < Xh[2]:
                Xs = (S(X0[0], Xc[0]), S(X0[1], Xc[1]), Z(S(X0[0], Xc[0]), S(X0[1], Xc[1])))
            else:
                Xs = (S(Xh[0], Xc[0]), S(Xh[1], Xc[1]), Z(S(Xh[0], Xc[0]), S(Xh[1], Xc[1])))
           #bpy.ops.mesh.primitive_cylinder_add(vertices=5, radius=0.002, depth=0.1, location=(Xs))
            if Xs[2] < Xh[2]:
                Xh = Xs
            else:
                Xh = (half(Xh[0], Xl[0]), half(Xh[1], Xl[1]), Z(half(Xh[0], Xl[0]), half(Xh[1], Xl[1])))
                Xg = (half(Xg[0], Xl[0]), half(Xg[1], Xl[1]), Z(half(Xg[0], Xl[0]), half(Xg[1], Xl[1])))
        else:
            break
        # bpy.ops.mesh.primitive_cylinder_add(vertices=3, radius=0.2, depth=0.1, location=(Xh))

        V = ((Xh[0], Xh[1], Xh[2]), (Xg[0], Xg[1], Xg[2]), (Xl[0], Xl[1], Xl[2]))
        midle = (Xh[2] + Xg[2] + Xl[2]) / 3
        Q = sqrt(pow((Xh[2] - midle), 2) + pow((Xg[2] - midle), 2) + pow((Xl[2] - midle), 2))
        N += 1
    #bpy.ops.mesh.primitive_cylinder_add(vertices=4, radius=0.005, depth=0.1, location=(Xl))

    face = ()  # кортеж индексов вершин одной грани
    FACES = ()  # кортеж кортежей индексов всех граней
    i = 0
    # Обход списка вершин

    while i < len(VERTEXES):
        face = (i, i + 1, i + 2)
        FACES += face,
        i += 3  # инкремент индекса, кратный длине строки
    # ----------------------------------------------------------------------------
    # Создание mesh-объекта
    #
    # создаём полигональную решётку
    # mesh = bpy.data.meshes.new(name='NEW_MESH_OBJECT')
    # # создаём объект, связанный с созданной полигональной решёткой
    # mesh_obj = bpy.data.objects.new('NEW_MESH_OBJECT', mesh)
    # # связываем объект с текущей сценой, полученной из контекста
    # scn = bpy.context.scene
    # scn.objects.link(mesh_obj)
    # scn.objects.active = mesh_obj
    # mesh_obj.select = True  # выделение объекта после его создания
    # # Наполняем mesh-объект ранее сгенерированной геометрией
    # mesh.from_pydata(VERTEXES, [], FACES)
    # # обновляем mesh-объект согласно его новым данным
    # mesh.update()
    print(VERTEXES[(len(VERTEXES)-4):(len(VERTEXES)-1)])
    print(" 1 метод завершил выполнение")


def f2():
    print("Вы во втором методе.  Результаты:")
    def Line(H, s):
        V = ()
        R = ()
        for i in range(0, len(H)):
            V += (t0 + i * dt, H[i], 10),
        for i in range(0, len(V) - 1):
            R += (i, i + 1),
        # mesh = bpy.data.meshes.new(name=s)
        # mesh_obj = bpy.data.objects.new((s), mesh)
        # scn = bpy.context.scene
        # scn.objects.link(mesh_obj)
        # scn.objects.active = mesh_obj
        # mesh.from_pydata(V, R, [])
        # mesh.update()

    def Zi(k, U, V):
        return Z[k - 1] + dt * (U + dt * V)

    def Xi(k, V):
        return X[k - 1] + dt * (Zi(k, u0, V))

    def setZ(U):
        for i in range(1, N + 1):
            Z[i] = (Zi(i, U, v0))

    def setX(V):
        for i in range(1, N + 1):
            X[i] = (Xi(i, V))

    x0 = 0
    xt = 1
    t0 = 0
    tt = 1

    z0 = 0
    zt = 1
    u0 = 1
    v0 = 1
    ez = 0.00001
    ex = 0.00001

    h = 0.001
    t = t0
    N = 10
    dt = (tt + t0) / (N)
    Z = [z0]
    X = [x0]
    for i in range(1, N + 1):
        Z.append(0)
        X.append(0)

    setZ(u0)
    z = Z[N]

    if abs(z - zt) > ez:

        setZ(u0 + h)
        zp = Z[N]
        setZ(u0 - h)
        zm = Z[N]

        if abs(zp - zt) < abs(zm - zt):
            print("Xeq")
            u0 += h
            setZ(u0)
            while abs(Z[N] - zt) > ez:
                u0 += h
                setZ(u0)
            u0 -= h
            setZ(u0)
        elif abs(zp - zt) > abs(zm - zt):
            u0 -= h
            setZ(u0)
            while abs(Z[N] - zt) > ez:
                u0 -= h
                setZ(u0)
            u0 += h
            setZ(u0)

    setX(v0)
    x = X[N]

    if abs(x - xt) > ex:
        setX(v0 + h)
        xp = X[N]
        setX(v0 - h)
        xm = X[N]
        if abs(xp - xt) < abs(xm - xt):
            v0 += h
            setX(v0)
            while abs(X[N] - xt) > ex:
                v0 += h
                setX(v0)
            v0 -= h
            setX(v0)
        elif abs(xp - xt) > abs(xm - xt):
            v0 -= h
            setX(v0)
            while abs(X[N] - xt) > ex:
                v0 -= h
                setX(v0)
            v0 += h
            setX(v0)
    print("U = " + str(u0) + "; V = " + str(v0) + "  : ")
    print('Z=',Z)
    print('X=',X)

    Line(Z, "z")
    Line(X, "X")
    print(" 1 метод завершил выполнение")

allMetod =['Cимплексный метод', 'Метод пристрелки']
methods = {1: 'Cимплексный метод',2: 'Метод пристрелки'}
while True :
    print('Доступные методы:', end=' ')
    for element in allMetod:
        print(element, end=', ')
    print()
    a=(int)(input())
    print(methods.get(a), end=":  ")
    if a==1:
        f()
    elif a==2:
        f2()
    else:
        print('Нет метода')
    print()

