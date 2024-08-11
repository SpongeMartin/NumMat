using Main.DN1

H = Hessenberg([2 7 2;
                3 2 7;
                0 1 4])

Q,R = qr(H)
Q
R
newH = qrIteration(H,1500)
newH