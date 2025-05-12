%Fun??es e f?rmulas do Poder Domínio NOMA
%Implementação das ferramentas para calcular o Sum Rate NOMA

%==================== BIBLIOTECAS =============================
pkg load statistics   #Versão 6.1.0
%==================== TABELA DE PARÂMETROS ====================
Phi = 70;        #graus -> ângulo de radiação (?)
T_csi = 1.0;     #Ganho do filtro óptico
f = 1.5;         #?nice refrativo
Csi_fov = 85;    #Graus -> ângulo fixado
d = 2.5;         #metros -> distância de linkagem
Apd = 1.0;       #centímetro quadrado -> área do fotodetector
ro_RIS = 0.95;   #reflectividade do RIS
%==============================================================
%=================== Parâmetros de TESTE ======================
#Utilizados no cálculo do cosseno (somente para ver se funciona)
x1 = 0;   #acess point
y1 = 12;  #acess point
z1 = 1;   #acess point
x2 = 2;   #user
y2 = 1;   #user
z2 = 0;   #user
%==============================================================
%============= Parâmetros que Precisam ser Definidos ==========
Ak = 1;     #área do k-ésimo espelho
An = 1;             #Área do n-ésimo elemento refrativo - to considerando o mesmo tamanho do espelho
#cosseno_phi_an = ; #trocar os parâmetros do cosseno_csi_KU
#cosseno_csi_an = ; #mesma ideia do anterior
#cosseno_phi_nu = ; #escrever expressão
#cosseno_csi_nu = ; #trocar parâmetros
%==============================================================
%=================== Definições ===============================
#Colocar as definições dos ângulos aqui
yaw = 45      #Extraído do artigo Modeling the Random Orientation of Mobile Devices
roll = 10     #Extraído do artigo Modeling the Random Orientation of Mobile Devices

#-----------------------------------------------------------------------------------
#alpha = -180 + (180 - (-180))*rand()#Distribuição Uniforme -> Ângulo Polar
#beta = laplace_rnd(0, 1)            #Média = 0//Parâmetro da Escala = 1
#while !(beta >= 0 && beta <= 90)    #Enquanto for menor < 0 ou > 90 gera novamente
#  beta = laplace_rnd(0, 1)          #Distribuição Laplace -> Ângulo Azimutal
#endwhile
#Aqui eu gerei os ângulos com as fórmulas acima
#e travei eles para variar os outros parâmetros
#(posição do ponto de acesso e do usuário)
alpha = 80.890
beta = 2.3823
#-----------------------------------------------------------------------------------
#Usuários - Coordenada z é sempre zero
u1 = [2, 4, 0];
u2 = [5, 2, 0];
u3 = [7.98, 4, 0];
u4 = [6.5, 8, 0];
u5 = [5.5, 4, 0];
u6 = [-2, 1, 0];
u7 = [-3, 4, 0];
u8 = [-3, 1, 0];
UxRoom1 = [2, 5, 7.98, 6.5, 5.5];
UyRoom1 = [4, 2, 4, 8, 4];
UxRoom2 = [-2, -3, -3];
UyRoom2 = [1, 4, 1];
Usuarios = {
      u1,
      u2,
      u3,
      u4,
      u5,
      u6,
      u7,
      u8
};
#-----------------------------------------------------------------------------------
VLC = [5, 5, 3];             #Coordenadas do vlc
MirrorArray = [0, 0, 2];     #Coordenadas do vetor de espelhos
#-----------------------------------------------------------------------------------
cossenos_csi = [];
cossenos_phi_KU = [];         #ângulo de irradiância para o caminho que sai da superfície k até o usuário u
cossenos_csi_AK = [];         #ângulo incidente para o caminho que sai do AP até a k-ésima superfície refletora
cossenos_csi_KU = [];         #ângulo incidente para o caminho que sai do k-ésimo espelho até o usuário u
cosseno_phi_AK = 0;           #ângulo de radiação para o caminho do AP até a k-ésima superfície refletora
#OBS: Como estamos considerando apenas um espelho, entendi que esse cosseno não muda (cálculo l:136)
#-----------------------------------------------------------------------------------
dKU = [];                                #distância do k-ésimo espelho até o usuario RIS u
d_UK = [];                               #distância entre o usuário e o espelho (superfície refletora) - room2
dAK = sqrt(sum((VLC - MirrorArray).^2)); #distância do AP até o k-ésimo espelho
#-----------------------------------------------------------------------------------
#Ganhos de canal
HLos = [];          #Ganho do Canal para o Caminho da Linha de Visão
Hnlos = [];         #Ganho do Canal sem a linha de visão
%==============================================================

#Ordem Lambertiana de emissão (m)
#Phi denota o ângulo de radiação
m = -(log2(cos(Phi)))^(-1)

#Ganho Óptico do concentrador não imaginário
G_csi = f^2/sin(Csi_fov)^2

#Distância entre o espelho e cada usuário u (room1)
for i = 1:5
  dKU(i) = sqrt(sum((Usuarios{i} - MirrorArray).^2));
endfor
dKU

#Distância entre a surpefície refletora e o usuário (room2) - como no artigo está separado,
#achei melhor separar aqui também
for i = 6:8
  d_UK(i) = sqrt(sum((Usuarios{i} - MirrorArray).^2));
endfor

#Cálculo do cosseno do ângulo de incidência
#As coordenadas não foram explicitadas no artigo - defini as coordenadas
# (x1, y1, z1) -> coordenadas do vetor posição do ponto de acesso (AP)
# (x2, y2, z2) -> coordenadas do vetor posição do usuário (u) -> z2 é sempre zero aqui
function Cosseno_csi = cosseno_csi(x1, y1, z1, x2, y2, z2, alpha, beta, d)
Cosseno_csi = ((abs(x1 - x2)/d) * cosd(beta) * sind(alpha) +
              (abs(y1 - y2)/d) * sind(beta) * cosd(alpha) +
              (abs(z1 - z2)/d) * cosd(alpha));
endfunction

#Cálculo de cosseno_phi
function Cosseno_phi = cosseno_phi(x_espelho, y_espelho, z_espelho, x_user, y_user, z_user, yaw, roll, d_UK)
Cosseno_phi = ((abs(x_espelho - x_user)/d_UK)*sind(yaw)*cosd(roll) +
               (abs(y_espelho - y_user)/d_UK)*cosd(yaw)*cosd(roll) +
               (abs(z_espelho - z_user)/d_UK)*sind(roll));
endfunction

#calcula o cosseno_csi e o cosseno_phi para cada usuário u no ->Room1<-
for i = 1:5
  cossenos_csi(i) = cosseno_csi(VLC(1), VLC(2), VLC(3), UxRoom1(i), UyRoom1(i), 0, alpha, beta, d);
  cossenos_csi_AK(i) = cosseno_csi(VLC(1), VLC(2), VLC(3), MirrorArray(1), MirrorArray(2), MirrorArray(3), alpha, beta, d);
  cossenos_csi_KU(i) = cosseno_csi(MirrorArray(1), MirrorArray(2), MirrorArray(3), UxRoom1(i), UyRoom1(i), 0, alpha, beta, d);
  cossenos_phi_KU(i) = cosseno_phi(MirrorArray(1), MirrorArray(2), MirrorArray(3), UxRoom1(i), UyRoom1(i), 0, yaw, roll, dKU(i));
endfor

cossenos_csi
cossenos_csi_AK
cossenos_csi_KU
cossenos_phi_KU

#Calcula cosseno_phi_AK (valor único na nossa situação)
cosseno_phi_AK = cosseno_phi(VLC(1), VLC(2), VLC(3), MirrorArray(1), MirrorArray(2), MirrorArray(3), yaw, roll, dAK)


#Ganho do Canal para o Caminho da Linha de Visão (HLos)
#Consideramos que o ângulo de incidência csi está dentro do intervalo correto
# 0 <= Csi <= Csi_fov
#Caso não esteja, Hlos = 0
function DevolveHlos = Hlos(m, G_csi, T_csi, Phi, cosseno_csi, Apd, d)
  DevolveHlos = ((((m + 1) * Apd)/((2 * pi) * (d^2))) * G_csi * T_csi * (cos(Phi)^m) * cosseno_csi);
endfunction

#Calcula HLos para cada usuário (o que varia é o cosseno_csi)
for i = 1:5
  HLos(i) = Hlos(m, G_csi, Phi, T_csi, cossenos_csi(i), Apd, d);
endfor

HLos
%-------------------------------------------------------------------------------------------------------------------------------------------

#Ganho do Canal para o Caminho SEM Linha de Visão (HNlos)
#Com superfícies refletoras
#Consideramos que o ângulo de incidência csi está dentro do intervalo correto
#G_csi e T_csi já estão definidos
function DevolveHNlos = HNlos(m, ro_RIS, Apd, dAK, dKU, Ak, G_csi, T_csi, cosseno_phi_AK, cosseno_csi_AK, cosseno_phi_KU, cosseno_csi_KU)
  DevolveHNlos = (ro_RIS * (((m + 1) * Apd)/(2*(pi^2) * dAK^2 * dKU^2)) * Ak * G_csi * T_csi *
                  cosseno_phi_AK^m * cosseno_csi_AK * cosseno_phi_KU * cosseno_csi_KU)
endfunction

#HNlos para cada usuário - Room1
for i = 1:5
  Hnlos(i) = HNlos(m, ro_RIS, Apd, dAK, dKU(i), Ak, G_csi, T_csi, cosseno_phi_AK, cossenos_csi_AK(i), cossenos_phi_KU(i), cossenos_csi_KU(i));
endfor
Hnlos
%-------------------------------------------------------------------------------------------------------------------------------------------

#Ganho do Canal SEM linha de visão (HNlosn)
#Ganho de canal para o para o caminho NLos para a propagação do sinal do AP
#através da n-ésima superfície REFRATORA e o u-ésimo usuário



%-------------------------------------------------------------------------------------------------------------------------------------------




