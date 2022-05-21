#Implementation of the Single-Cut L-Shaped Algorithm for Finite Support R.V.
# LINMA2491 - Operations Research - Homework 2
# Group Members: Phong Do, Julien Donadello, Matthias Geurts, Eduardo Vannini
# Date: 20-04-2021

########################### Parameters ##############################

# Number of the test, choose among 1, 2 and 3
global TEST_NUMBER = 4

######################################################################

include("lshape.jl")


# Problem parameters

function Test_1()
    c = [100; 150]
    p = [0.4; 0.6]
    q = [-24 -28; -28 -32]
    K = 2

    A_eq = [0 0]
    b_eq = [0]
    A_ge = [1 0; 0 1]
    b_ge = [40; 20]
    A_le = [1 1]
    b_le = [120]

    W_eq = [0 0]
    W_ge = [1 0; 0 1]
    W_le = [6 10; 8 5; 1 0; 0 1]

    T_eq = zeros(K, 1, 2)
    T_eq[1, :, :] = [0 0]
    T_eq[2, :, :] = [0 0]

    T_ge = zeros(K, 2, 2)
    T_ge[1, :, :] = [0 0; 0 0]
    T_ge[2, :, :] = [0 0; 0 0]

    T_le = zeros(K, 4, 2)
    T_le[1, :, :] = [-60 0; 0 -80; 0 0; 0 0]
    T_le[2, :, :] = [-60 0; 0 -80; 0 0; 0 0]

    h_eq = [0; 0]
    h_ge = [0 0; 0 0]
    h_le = [0 0 500 100; 0 0 300 300]

    return c,p,q,K,A_eq,b_eq,A_ge,b_ge,b_le,W_eq,W_ge,W_le,A_le,T_eq,T_ge,T_le,h_eq,h_ge,h_le
end

function Test_2()
    c = [0]
    p = [1/3; 1/3; 1/3]
    q = [1 1; 1 1; 1 1]
    K = 3

    A_eq = reshape([0], 1, 1)
    b_eq = [0]
    A_ge = reshape([1], 1, 1)
    b_ge = [0]
    A_le = reshape([1], 1, 1)
    b_le = [10]

    W_eq = [1 -1]
    W_ge = [1 0; 0 1]
    W_le = [0 0]

    T_eq = zeros(K, 1, 1)
    T_eq[1, :, :] = [1]
    T_eq[2, :, :] = [1]
    T_eq[3, :, :] = [1]

    T_ge = zeros(K, 2, 1)
    T_ge[1, :, :] = [0; 0]
    T_ge[2, :, :] = [0; 0]
    T_ge[3, :, :] = [0; 0]

    T_le = zeros(K, 1, 1)
    T_le[1, :, :] = [0]
    T_le[2, :, :] = [0]
    T_le[3, :, :] = [0]

    h_eq = [1; 2; 4]
    h_ge = [0 0; 0 0; 0 0]
    h_le = [0; 0; 0]

    return c,p,q,K,A_eq,b_eq,A_ge,b_ge,b_le,W_eq,W_ge,W_le,A_le,T_eq,T_ge,T_le,h_eq,h_ge,h_le
end

function Test_3()
    c = [3; 2]
    p = [1/4; 1/4; 1/4; 1/4]
    q = [-15 -12; -15 -12; -15 -12; -15 -12]
    K = 4

    A_eq = [0 0]
    b_eq = [0]
    A_ge = [1 0; 0 1]
    b_ge = [0; 0]
    A_le = [0 0]
    b_le = [0]

    W_eq = [0 0]
    W_ge = [1 0; 0 1; 1 0; 0 1]
    W_le = [3 2; 2 5; 1 0; 0 1]

    T_eq = zeros(K, 1, 2)
    T_eq[1, :, :] = [0 0]
    T_eq[2, :, :] = [0 0]
    T_eq[3, :, :] = [0 0]
    T_eq[4, :, :] = [0 0]

    T_ge = zeros(K, 4, 2)
    T_ge[1, :, :] = [0 0; 0 0; 0 0; 0 0]
    T_ge[2, :, :] = [0 0; 0 0; 0 0; 0 0]
    T_ge[3, :, :] = [0 0; 0 0; 0 0; 0 0]
    T_ge[4, :, :] = [0 0; 0 0; 0 0; 0 0]

    T_le = zeros(K, 4, 2)
    T_le[1, :, :] = [-1 0; 0 -1; 0 0; 0 0]
    T_le[2, :, :] = [-1 0; 0 -1; 0 0; 0 0]
    T_le[3, :, :] = [-1 0; 0 -1; 0 0; 0 0]
    T_le[4, :, :] = [-1 0; 0 -1; 0 0; 0 0]

    h_eq = [0; 0; 0; 0]
    h_ge = [4.8 6.4 0 0; 3.2 6.4 0 0; 3.2 3.2 0 0; 4.8 3.2 0 0]
    h_le = [0 0 6 8; 0 0 4 8; 0 0 4 4; 0 0 6 4]

    return c,p,q,K,A_eq,b_eq,A_ge,b_ge,b_le,W_eq,W_ge,W_le,A_le,T_eq,T_ge,T_le,h_eq,h_ge,h_le
end
function Test_4()
    c = [150; 230; 260]
    p = [1/3; 1/3; 1/3]
    q = [238 210 -170 -150 -36 -10; 238 210 -170 -150 -36 -10; 238 210 -170 -150 -36 -10]
    K = 3

    A_eq = [0 0 0]
    b_eq = [0]
    A_ge = [1 0 0; 0 1 0; 0 0 1]
    b_ge = [0.; 0.; 0.]
    A_le = [1 1 1]
    b_le = [500]

    W_eq = [0 0 0 0 0 0]
    T_eq = zeros(K, 1, 3)
    T_eq[1, :, :] = [0 0 0]
    T_eq[2, :, :] = [0 0 0]
    T_eq[3, :, :] = [0 0 0]
    h_eq = [0; 0; 0]

    W_ge = [1 0 -1 0 0 0; 0 1 0 -1 0 0]
    # W_ge[1, :, :] = [1 0 -1 0 0 0; 0 1 0 -1 0 0]
    # W_ge[2, :, :] = [1 0 -1 0 0 0; 0 1 0 -1 0 0]
    # W_ge[3, :, :] = [1 0 -1 0 0 0; 0 1 0 -1 0 0]
    T_ge = zeros(K, 2, 3)
    T_ge[1, :, :] = [3.0 0 0; 0 3.6 0]
    T_ge[2, :, :] = [2.5 0 0; 0 3.0 0]
    T_ge[3, :, :] = [2.0 0 0; 0 2.4 0]
    h_ge = [200 240; 200 240; 200 240]

    W_le = [0 0 0 0 1 1; 0 0 0 0 1 0] #zeros(K,2,6) 
    # W_le[1, :, :] = [0 0 0 0 1 1; 0 0 0 0 1 0]
    # W_le[2, :, :] = [0 0 0 0 1 1; 0 0 0 0 1 0]
    # W_le[3, :, :] = [0 0 0 0 1 1; 0 0 0 0 1 0]
    T_le = zeros(K, 2, 3)
    T_le[1, :, :] = [0 0 -24; 0 0 0]
    T_le[2, :, :] = [0 0 -20; 0 0 0]
    T_le[3, :, :] = [0 0 -16; 0 0 0]
    h_le = [0 6000; 0 6000; 0 6000]

    return c,p,q,K,A_eq,b_eq,A_ge,b_ge,b_le,W_eq,W_ge,W_le,A_le,T_eq,T_ge,T_le,h_eq,h_ge,h_le
end

if TEST_NUMBER == 1
    c,p,q,K,A_eq,b_eq,A_ge,b_ge,b_le,W_eq,W_ge,W_le,A_le,T_eq,T_ge,T_le,h_eq,h_ge,h_le = Test_1()
elseif TEST_NUMBER == 2
    c,p,q,K,A_eq,b_eq,A_ge,b_ge,b_le,W_eq,W_ge,W_le,A_le,T_eq,T_ge,T_le,h_eq,h_ge,h_le = Test_2()
elseif TEST_NUMBER == 3
    c,p,q,K,A_eq,b_eq,A_ge,b_ge,b_le,W_eq,W_ge,W_le,A_le,T_eq,T_ge,T_le,h_eq,h_ge,h_le = Test_3()
else
    c,p,q,K,A_eq,b_eq,A_ge,b_ge,b_le,W_eq,W_ge,W_le,A_le,T_eq,T_ge,T_le,h_eq,h_ge,h_le = Test_4()
end

x_opt = L_Shape_Algorithm(c,p,q,K,A_eq,b_eq,A_ge,b_ge,b_le,W_eq,W_ge,W_le,A_le,T_eq,T_ge,T_le,h_eq,h_ge,h_le)
@printf("Optimal solution is x_opt = %s\n", x_opt)