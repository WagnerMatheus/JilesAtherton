%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Implementação do modelo de Jiles-Atherton escalar [1] (e, por
% consequência isotrópico) para cálculo da histerese em materiais
% magnéticos.
%
% EEL510236-41000056DO/ME (20251) - Materiais Elétricos e Magnéticos
%
% Matheus Henrique Wagner
%
%   De acordo com o modelo de Jiles-Atherton, a histerese magnética é modelada
% através de uma Equação Diferencial Ordinária (EDO) no tempo (t). A variável a
% ser integrada (ou variável dinâmica, também às vezes chamada de variável de
% estado) é a magnetização irreversível (M_ir) do material. O modelo também
% pode ser formulado para que a magnetização total (M) seja a variável de
% estado.
%
%   Não houve problemas com essa formulação no problema escalar, ao menos nos
% casos testados com magnetização inicial nula. Já no caso vetorial
% (Bergqvist) [1], algumas questões numéricas fazem com que a formulação
% baseada em `M_ir` seja cerca de 500x mais eficiente.
%
%   Além disso, levando em consideração que em modelos de elementos finitos a
% variável de estado é comumente o potencial (vetor) magnético, de onde se
% encontra de forma direta a indução magnética, a formulação aqui utilizada
% toma como excitação a indução, diferentemente do modelo original, que é
% excitado por campo magnético.
%
%   OBS.: este script foi testado em Octave e pode necessitar pequenas
% adaptações de sintaxe para rodar em MATLAB. Uma alternativa é utilizar o
% `OctaveOnline`:
% https://octave-online.net/workspace~afrqkndQgtEXNgCGTdPQCnxbcHehzbjhXwikwxfcRFZpyxVj
%
%   1. FORMULAÇÃO
%
%   Neste script a formulação com `M` sendo a variável de estado é utilizada.
% Ou seja, o problema a ser resolvido é
%
%                          - t
%                        /
%                       |  dM
%     M(t) = M(ti) +    |  -- dt ,  (1)
%                       |  dt
%                      /
%                    -  ti
%
% sendo `ti` o tempo inicial, e utiliza-se a indução magnética (B) como
% excitação
%
%     dM   dM dB
%     -- = -- -- ,  (2)
%     dt   dB dt
%
%               dM
%               --
%     dM        dH
%     -- = -------------- ,
%     dB             dM
%          mu0 ( 1 + -- )
%                    dH
%
%                   dM_an         dM_ir
%                 c ----- + (1-c) -----
%     dM             dHe           dHe
%     -- = ------------------------------------ ,
%     dH                dM_an          dM_ir
%          1 - alpha (c ------ + (1-c) ----- )
%                        dHe            dHe
%
% sendo `alpha` (acoplamento inter-domínio) e `c` (reversibilidade da
% magnetização) parâmetros do modelo, `mu0` a permeabilidade magnética do
% vácuo, `H` o campo magnético e (baseando-se em Bergqvist [1] para a
% derivada de `M_ir`)
%
%     He = H + alpha M,
%
%     M_an = ms (coth(He/a)-(a/He)),
%
%                 -  | M_an - M_ir |                    dH
%     dM_ir     /    --------------- , se (M_an - M_ir) -- > 0     (3)
%     ----- = -             k                           dt      ,
%      dHe      \
%                 -         0        , caso contrário
%
% sendo `ms` (magnetização de saturação), `a` (densidade de paredes de domínio)
% e `k` (energia média necessária para romper `fixamento` (pinning)) são
% parâmetros do modelo e
%
%            M - c M_an
%     M_ir = ---------- .  (4)
%              1-c
%
%   Neste trabalho, `dH/dt` na equação 3 será substituído por `dB/dt`, visto
% que a derivada de ambos no caso escalar tem o mesmo sinal, assim
%
%                 -  | M_an - M_ir |                    dB
%     dM_ir     /    --------------- , se (M_an - M_ir) -- > 0
%     ----- = -             k                           dt      .
%      dHe      \
%                 -         0        , caso contrário
%
%   Também nota-se que na equação 4 fica implícita a relação entre as
% componentes reversível (M_an) e irreversível (M_ir) da magnetização total
% (M):
%
%     M = c M_an + (1-c) M_ir.
%
%   OBS.: os parâmetros do modelo de histerese utilizados para simulação neste
% script provém de [2], tabela 6.1, coluna `Conhecido`. São eles:
%
%   +----------+-------------+
%   | ms [A/m] | 1,47 x 10^6 |
%   +----------+-------------+
%   | k [A/m]  | 70,0        |
%   +----------+-------------+
%   | c        | 340 x 10^-3 |
%   +----------+-------------+
%   | a [A/m]  | 89,0        |
%   +----------+-------------+
%   | alpha    | 169 x 10^-6 |
%   +----------+-------------+
%
%   2. INTEGRAÇÃO
%
%   A resolução numérica de EDOs no tempo (também chamadas de Problema de
% Valor Inicial (PVI)) é realizada através de integradores.
%
%   Esse tópico é tema de pesquisa há várias décadas, pois é fundamental para
% muitas áreas da física e engenharia. Na elétrica, qualquer estudo transiente
% de equipamento ou circuito precisa resolver numericamente um PVI. Dessa
% forma, existem diversos algoritmos exaustivamente estudados e testados à
% nossa disposição.
%
%   Para este script em questão, é relevante saber que os algoritmos
% implementados para resolução de PVIs ditos explícitos (caso do modelo de
% Jiles-Atherton) são padronizados para ter como `entrada` a função f(x, t)
% definida por
%
%               dx
%     f(x, t) = -- ,
%               dt
%
% sendo `x` a variável de estado, e também o valor inicial (x_0) e o intervalo
% de integração [ti, tf].
%
%   3. IMPLEMENTAÇÃO
%
%   A implementação consiste em fornecer uma função `dmdt` que calcule o lado
% esquerdo da equação 2, tendo como entradas o valor atual de `M` e a excitação
% `B` e `dB/dt`, ou seja,
%
%     f(M, t) = dmdt(M, B(t), dBdt(t)).
%
%   Assim, também se faz necessário fornecer uma excitação de indução
%   magnética arbitrária, juntamente com uma função para cálculo da sua
%   derivada.
%
%   3.1. INTEGRADORES DE EDO DO OCTAVE/MATLAB
%
%   Runge-Kutta
%   -----------
%
%   Os métodos Runge-Kutta são implementados pela própria biblioteca do Octave.
%
%     `ode45`: integra um sistema de EDOs non-stiff (ou EDAs (Equações
%     Diferencias Algébricas) de indíce 1) usando o método de Dormand-Prince de
%     ordem alta e passo variável
%
%     `ode23`: integra um sistema de EDOs non-stiff (ou EDAs de índice 1).
%     Utiliza o método de Bogacki-Shampine de ordem 3 e adapta o tamanho de
%     passo local de acordo com uma tolerância especificada pelo usuário.
%
%     `ode23s`: integra um sistem de EDOs stiff (ou EDAs de índice 1)
%     utilizando um método de Rosenbrock de ordem 2 modificado.
%
%
%   Métodos de passo múltiplo
%   -------------------------
%
%   Métodos que utilizam o IDA, da biblioteca (c) SUNDIALS:
%
%     `ode15s`: integra um sistema de EDOs stiff (ou EDAs de índice 1)
%     utilizando um método de passo e ordem variáveis baseado nas fórmulas de
%     diferenciação atrasadas (Backward Difference Formulas (BDF)).
%
%     `ode15i`: integra um sistema de EDOs ímplicitas (ou EDAs de índice 1)
%     utilizando o mesmo método da `ode15s`.
%
%   Método que utiliza o LSODE, da biblioteca (Fortran 77) ODEPACK:
%
%     `lsode`: integra um sistema de EDOs stiff ou non-stiff. Utiliza métodos
%     de passo e ordem variáveis, sendo (BDF) para EDOs stiff.
%
%   OBS.: neste script limita-se ao uso do solver `lsode`. A utilização de
% qualquer um dos outros solvers tem sintaxe praticamente idêntica.
%
%   4. VISUALIZAÇÃO
%
%   Alguns gráficos são gerados para visualização dos resultados, são eles:
%
%  - Laço de histerese;
%  - Indução no tempo;
%  - Campo magnético no tempo;
%  - Componentes da magnetização no tempo.
%
%   Esses gráficos são gerados para os seguintes casos:
%
%  - Curva de magnetização inicial;
%  - Aumento gradual da indução máxima alternada, formando laços de
%    variadas amplitudes;
%  - Indução com harmônicas de 3º e 5º grau, formando laços menores
%    deslocados do 0,0;
%
%   5. REFERÊNCIAS
%
% [1]: Gyselinck, Johan & Dular, Patrick & Sadowski, Nelson & Leite, J.V.
% & Bastos, J.P.A.. (2004). Incorporation of a Jiles-Atherton vector
% hysteresis model in 2D FE magnetic field computations: Application of
% the Newton-Raphson method. COMPEL: The International Journal for
% Computation and Mathematics in Electrical and Electronic Engineering.
% 23. 685-693. 10.1108/03321640410540601.
%
% [2]: Batistela, N. J. (2001). Caracterização e modelagem eletromagnética
% de lâminas de aço silício. Tese (doutorado) - Universidade Federal de
% Santa Catarina, Florianópolis, SC, Brasil.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

%%%%% CONSTANTES FÍSICAS %%%%%
global mu0 = 4e-7*pi;

%%%%% PARÂMETROS DO MATERIAL %%%%%
% parâmetros fornecidos em [2]
global ms = 1.47e6; global a = 89; global k = 70; global c = 340e-3;
global alpha = 169e-6;

% parâmetros fornecidos em [1]
%global ms = 1145500; global a = 59; global k = 99; global c = .55;
%global alpha = 1.3e-4;

% fronteira do intervalo [-He/a = -man_iso_lc, He/a = man_iso_lc] dentro do
% qual M_an é linearmente continuada, evitando divisão por números próximos de
% zero
global man_iso_lc = 1e-4;
global slope_iso_lc = ms * (1.0/tanh(man_iso_lc)-(1/man_iso_lc))/(a*man_iso_lc);

%%%%% DEFINIÇÃO DAS FUNÇÕES DO MODELO %%%%%
% M_an
function result = man_iso(He)
  global a; global man_iso_lc; global ms; global slope_iso_lc;
  if abs(He/a) > man_iso_lc
    % região longe de zero
    result = ms * (1.0/tanh(He/a)-(a/He));
  else 
    % região próxima de zero
    result = slope_iso_lc * He;
  end
end

% dM_an
% -----
%  dHe
function result = dmandHe(He)
  global a; global man_iso_lc; global ms; global slope_iso_lc;
  if abs(He/a) > man_iso_lc
    % região longe de zero
    result = ms * (a/(He*He) - (1.0/sinh(He/a))^2 / a);
  else
    % região próxima de zero
    result = slope_iso_lc;
  end
end

% dM_ir
% -----
%  dHe
% definição baseada em Bergqvist [1]
function result = dmirdHe(man, mi, dBdt)
  global k;
  aux = (man - mi)/k;
  if dBdt * aux > 0
    result = abs(aux);
  else
    result = 0;
  end
end

% dM
% --
% dt
function result = dmdt(M, B, dBdt)
  global mu0; global alpha; global c;
  H = B/mu0 - M;
  He = H + alpha * M;
  M_an = man_iso(He);
  M_ir = (M - c*M_an) / (1-c);
  dM_andHe = dmandHe(He);
  dM_irdHe = dmirdHe(M_an, M_ir, dBdt);
  dMdH = (c*dM_andHe + (1-c)*dM_irdHe) ...
          / (1 - alpha * (c*dM_andHe + (1-c)*dM_irdHe));
  dMdB = dMdH / (mu0 * (1 + dMdH));
  result = dMdB * dBdt;
end

%%%%% DEFINIÇÃO DAS FUNÇÕES DE EXCITAÇÃO %%%%%

% intervalo de simulação
global ti = 0; global tf = 1;

% indução máxima (para algumas simulações)
global B_max = 1.5;

% B(t)
function result = B(t)
  % São disponibilizadas diversas curvas `B(t)` para geração de cada um
  % dos gráficos apresentados no trabalho. Para simular, basta descomentar
  % a curva desejada e comentar todas as outras. Não esquecer de fazer o
  % mesmo para dBdt(t).

  global ti; global tf; global B_max;

  % Curva de magnetização inicial
  %result = t * B_max/(tf-ti);

  % `n-1` laços menores (e 1 externo) (plotar com no mínimo 100n pontos)
  %n = 8; T = (tf-ti)/n;
  %j = floor((t-ti)/T) + 1;
  %result = j/n * 1.4*sin(2*pi*(2/T)*t);

  % harmônica de nº grau com amplitude relativa h_amp e defasagem phi
  % caso n = 3
  n = 3; h_amp = .5; phi = 30 * pi/180;
  result = B_max/1.3 * (sin(2*pi*(2/(tf-ti))*t) ...
    + h_amp * sin(2*pi*(n*2/(tf-ti))*t + phi));

  % caso n = 5
  %n = 5; h_amp = .5; phi = 80 * pi/180;
  %result = 1/1.6 * B_max*(sin(2*pi*(2/(tf-ti))*t) ...
  %  + h_amp * sin(2*pi*(n*2/(tf-ti))*t + phi));
end

% dB
% --
% dt
function result = dBdt(t)
  global ti; global tf; global B_max;

  % magnetização inicial
  %result = B_max/(tf-ti);

  % `n-1` laços menores (e 1 externo) (plotar com no mínimo 100n pontos)
  %n = 8; T = (tf-ti)/n;
  %j = floor((t-ti)/T) + 1;
  %result = j/n * 2*pi*(2/T) * 1.4*cos(2*pi*(2/T)*t);

  % harmônica de nº grau com amplitude relativa h_amp e defasagem phi
  % caso n = 3
  n = 3; h_amp = .5; phi = 30 * pi/180;
  result = B_max/1.3 * (2*pi*(2/(tf-ti)) * cos(2*pi*(2/(tf-ti))*t) ...
    + 2*pi*(n*2/(tf-ti)) * h_amp * cos(2*pi*(n*2/(tf-ti))*t + phi));

  % caso n = 5
  %n = 5; h_amp = .5; phi = 80 * pi/180;
  %result = 1/1.6 * B_max * (2*pi*(2/(tf-ti)) * cos(2*pi*(2/(tf-ti))*t) ...
  %  + 2*pi*(n*2/(tf-ti)) * h_amp * cos(2*pi*(n*2/(tf-ti))*t + phi));
end

%%%%% DEFINIÇÃO DOS PARÂMETROS E RESOLUÇÃO %%%%%

% magnetização inicial (tal que H(ti) = 0)
M_0 = B(ti)/mu0;

% número de passos no tempo para os quais se tem interesse na solução
num_output_steps = 300;

ode_options = odeset(
  'RelTol', 1e-4,
  'AbsTol', 1e-6,
  'MaxStep', 1e-1,
  'InitialStep', 1e-3);

f = @(M, t) dmdt(M, B(t), dBdt(t));

%%%%% SOLUÇÃO UTILIZANDO O LSODE %%%%%

sol_time = transpose(linspace(ti, tf, num_output_steps));
M_vec = lsode(f, M_0, sol_time);

%%%%% VISUALIZAÇÃO %%%%%

% cálculo dos vetores para plotagem
B_vec = zeros(size(sol_time));
H_vec = zeros(size(sol_time));
Man_vec = zeros(size(sol_time));
Mir_vec = zeros(size(sol_time));

for i=1:length(sol_time)
  B_vec(i) = B(sol_time(i));
  H_vec(i) = B_vec(i)/mu0-M_vec(i);
  He = H_vec(i) + alpha * M_vec(i);
  Man_vec(i) = man_iso(He);
  Mir_vec(i) = (M_vec(i) - c * Man_vec(i)) / (1-c);
end

% salvar valores em arquivo
dlmwrite("results.txt", [sol_time, H_vec, B_vec, Man_vec, Mir_vec], "delimiter", " ", "newline", "\n")

close all

set(0, "defaultlinelinewidth", 2, "defaultaxesfontsize", 12);

% plot do laço de histerese
figure
plot(H_vec, B_vec)
title("Laço de Histerese")
xlabel("H [A/m]")
ylabel("B [T]")

% plot do campo magnético no tempo
figure
plot(sol_time, H_vec)
title("Campo Magnético no Tempo")
xlabel("t [s]")
ylabel("H [A/m]")

% plot da indução no tempo
figure
plot(sol_time, B_vec)
title("Indução Magnética no Tempo")
xlabel("t [s]")
ylabel("B [T]")

% plot das componentes da magnetização no tempo
figure
hold on
plot(sol_time, Man_vec*1e-6)
plot(sol_time, Mir_vec*1e-6)
title("Componentes da Magnetização no Tempo")
xlabel("t [s]")
ylabel("[MA/m]")
legend("Mag. anisterética", "Mag. irrevesível")

