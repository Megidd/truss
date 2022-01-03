% This code is developed for analysis of 2D trusses
% Geometric Nonlinearity can be acounted but truss elements are linear
% elastic. Node cordinates have units of "m" and external forces should be
% defined in "N".
clc
clear
%% Truss Definition
% Nodes=input('input node cordinates in "m" with [x y] format.1st row belongs to the node 1:');
Nodes = [0	0
    3	0
    0	3
    3	3
    0	6
    3	6
    ];
% Elements=input('input element data in [Node-i Node-j A E] format. 1st row belongs to the element 1:');
Elements = [1	2	0.000302	2.00E+11
    1	3	0.000302	2.00E+11
    2	4	0.000302	2.00E+11
    3	2	0.000302	2.00E+11
    3	4	0.000302	2.00E+11
    3	5	0.000302	2.00E+11
    4	6	0.000302	2.00E+11
    3	6	0.000302	2.00E+11
    5	6	0.000302	2.00E+11

    ];

%Fext=input('input nodal external force in "N" with format of [Node Magnitude Agnle]:');
Fext = [3	50000	0
    5	111803.3989	-26.57
    6	50000	-90
    ];

%Supports=input('input Boundary conditions with [Node type Orientation]:');
Supports = [1	2	0
    2	2	0
    ];
%% Element Stiffness Calculation
for i = 1:size(Elements, 1)
    n1 = Elements (i, 1);
    n2 = Elements (i, 2);

    x1 = Nodes(n1, 1);
    x2 = Nodes(n2, 1);
    Dx = x2 - x1;

    y1 = Nodes(n1, 2);
    y2 = Nodes(n2, 2);
    Dy = y2 - y1;

    L(i) = (Dx^2 + Dy^2)^0.5; %Element Length
    Tet(i) = atan2(Dy, Dx); %Element Angle

    A = Elements(i, 3); %Cross-Section Area
    E = Elements(i, 4); %Elastic Modulus

    S = sin(Tet(i));
    C = cos(Tet(i));

    k{i} = E * A / L(i) * [C^2 S * C -C^2 -S * C;
                S * C S^2 -C * S -S^2;
                -C^2 -S * C C^2 S * C;
                -S * C -S^2 C * S S^2];
end

%% Global Stiffness Assemble
K = zeros(2 * size(Nodes, 1));
K0 = K;

for i = 1:size(Elements, 1)
    n1 = Elements(i, 1);
    n2 = Elements(i, 2);

    K0(2 * n1 - 1:2 * n1, 2 * n1 - 1:2 * n1) = k{i}(1:2, 1:2);
    K0(2 * n1 - 1:2 * n1, 2 * n2 - 1:2 * n2) = k{i}(1:2, 3:4);
    K0(2 * n2 - 1:2 * n2, 2 * n1 - 1:2 * n1) = k{i}(3:4, 1:2);
    K0(2 * n2 - 1:2 * n2, 2 * n2 - 1:2 * n2) = k{i}(3:4, 3:4);

    K = K + K0;
    K0(:, :) = 0;
end

%% External Forces

F0 = zeros(2 * size(Nodes, 1), 1);

for i = 1:size(Fext, 1)
    Fnode = Fext(i, 1);
    Fmag = Fext(i, 2);
    Fteta = Fext(i, 3) / 180 * pi;
    fx = Fmag * cos(Fteta);
    fy = Fmag * sin(Fteta);
    F0(2 * Fnode - 1:2 * Fnode, 1) = [fx; fy];
end

%% Boundary Conditions
q = 0;

for i = 1:size(Supports, 1);

    Snode = Supports(i, 1);
    Stype = Supports(i, 2);
    Sorien = Supports(i, 3);

    if Stype == 1 % if Roller

        if Sorien == 1 % if Horizontal
            q = q + 1;
            zeroDOF(q) = 2 * Snode - 1;
        elseif Sorien == 2 %if vertical
            q = q + 1;
            zeroDOF(q) = 2 * Snode;
        end

    elseif Stype == 2 %if Pin
        q = q + 2;
        zeroDOF(q - 1:q) = 2 * Snode - 1:2 * Snode; %list of DOFs with zero disp.
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Iteration start

%% Solving Equation

Kc = K;
Fc = F0;

Kc(:, zeroDOF) = []; %removing total stiffness columns with zero DOF
Kc(zeroDOF, :) = []; %removing total stiffness rows with zero DOF

Fc(zeroDOF, :) = []; %removing force rows with zero DOF

U0 = Kc^ - 1 * Fc; %Non-zero nodal displacements

u_all = 1:2 * size(Nodes, 1);
u_nonzero = u_all;
u_nonzero(zeroDOF) = [];

U(u_all, 1) = 0;
U(u_nonzero, 1) = U0; % Nodal displacement vector

F = K * U; % Nodal force vector

%% Element Forces and Elongations

for i = 1:size(Elements, 1);
    S = sin(Tet(i));
    C = cos(Tet(i));
    n1 = Elements(i, 1);
    n2 = Elements(i, 2);

    Delta = [-C, -S, C, S] * [U(2 * n1 - 1); U(2 * n1); U(2 * n2 - 1); U(2 * n2)];

    A = Elements(i, 3);
    E = Elements(i, 4);

    P = E * A / L(i) * Delta;
    ElResult(i, 1) = Delta; %Element elongation
    ElResult(i, 2) = P; %Element force
    ElResult(i, 3) = P / A; %Element stress
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Iteration end

%%Display Results

%fprintf('\n')
%for i=1:size(Elements,1);
%   fprintf('Element (%g) Elongation: %g mm\n',i,ElResult(i,1)*1000)
%end

fprintf('\n')

for i = 1:size(Elements, 1);
    fprintf('Element (%g) Force: %g kN\n', i, ElResult(i, 2))
end

fprintf('\n')
cntr = 1;

for i = 1:2:size(U, 1);
    fprintf('Node (%g) Displacement X: %g mm\n', cntr, U(i) * 1000)
    fprintf('Node (%g) Displacement Y: %g mm\n', cntr, U(i + 1) * 1000)
    cntr++;
end
