function [H, Gamma] = H(usuario, VLC, espelhos, refratores, yaw, roll, Eta_c, lambda, iota)
  %==================== BIBLIOTECAS =============================
  pkg load statistics   #Versão 6.1.0
  %==================== TABELA DE PARÂMETROS ====================
  num_esp = length(espelhos);           #Numero de espelhos
  num_ref = length(refratores);         #Numero de refratores
  ang_meia_pot = 70;                    #graus -> semi-ângulo de meia potencia
  T_csi = 1.0;                          #Ganho do filtro óptico
  f = 1.5;                              #?nice refrativo
  Csi_fov = 85;                         #Graus -> ângulo fixado
  d = 2.5;                              #metros -> distância de linkagem
  Apd = 1.0;                            #centímetro quadrado -> área do fotodetector
  ro_RIS = 0.95;                        #reflectividade do RIS
  m = -(log2(cosd(ang_meia_pot)))^(-1); #Ordem Lambertiana de emissão (m)
  G_csi = f^2/sind(Csi_fov)^2;          #Ganho Óptico do concentrador não imaginário
  Eta_a = 1.0;                          #indice refrativo do ar
  r_eff = 12;                           # pm/V -> electro-optic coefficient
  Vth = 1.34;                           # V -> critical voltage threshold
  Eta_o = 1.5;
  Eta_e = 1.7;
  D = 0.75;                             # mm -> profundidade da celula de cristal-liquido
  %==============================================================
  %============= Parâmetros que Precisam ser Definidos ==========
  Ak = 0.01;     #área do k-ésimo espelho
  An = 0.01;     #Área do n-ésimo elemento refrativo - to considerando o mesmo tamanho do espelho
  #-----------------------------------------------------------------------------------
  #alpha = -180 + (180 - (-180))*rand()#Distribuição Uniforme -> Ângulo Polar
  #beta = laplace_rnd(0, 1)            #Média = 0//Parâmetro da Escala = 1
  #while !(beta >= 0 && beta <= 90)    #Enquanto for menor < 0 ou > 90 gera novamente
  #  beta = laplace_rnd(0, 1)          #Distribuição Laplace -> Ângulo Azimutal
  #endwhile
  #Aqui eu gerei os ângulos com as fórmulas acima
  #e travei eles para variar os outros parâmetros
  #(posição do ponto de acesso e do usuário)
  alpha = 80.890;
  beta = 2.3823;
  #-----------------------------------------------------------------------------------
  cosseno_phi_AU = 0;        #ângulo de irradiância para o caminho que sai do AP até o usuário u
  cosseno_csi_AU = 0;        #ângulo incidente para o caminho que sai do AP até o usuário u

  cossenos_phi_KU = [];      #ângulo de irradiância para o caminho que sai da superfície k até o usuário u
  cossenos_csi_AK = [];      #ângulo incidente para o caminho que sai do AP até a k-ésima superfície refletora
  cossenos_csi_KU = [];      #ângulo incidente para o caminho que sai do k-ésimo espelho até o usuário u
  cossenos_phi_AK = [];      #ângulo de radiação para o caminho do AP até a k-ésima superfície refletora

  cossenos_phi_NU = [];      #cosseno do ângulo de irradiância para o caminho que sai da superfície n até o usuário u
  cossenos_csi_AN = [];      #cosseno do ângulo incidente para o caminho que sai do AP até a n-ésima superfície refratora
  cossenos_csi_NU = [];      #cosseno do ângulo incidente para o caminho que sai do n-ésimo espelho até o usuário u
  senos_csi_NU = [];         #seno do ângulo incidente para o caminho que sai do n-ésimo espelho até o usuário u
  cossenos_phi_AN = [];      #cosseno do ângulo de radiação para o caminho do AP até a n-ésima superfície refratora
  #teta eh o angulo depois de passar por refracao ao entrar no LC RIS
  senos_teta = [];
  cossenos_teta = [];
  #-----------------------------------------------------------------------------------
  dKU = [];                                #distância do k-ésimo espelho até o usuario u
  dNU = [];                                #distância do n-ésimo espelho até o usuario u
  dAK = []; #distância do AP até o k-ésimo espelho
  dAN = []; #distância do AP até o n-ésimo refrator (superfície refratora)
  #-----------------------------------------------------------------------------------
  #Ganhos de canal
  HLos = 0;           #Ganho do Canal para o Caminho da Linha de Visão
  Hnlosk = [];        #Ganho do Canal sem a linha de visão via superficies refletoras
  Hnlosn = [];        #Ganho do Canal sem a linha de visão via superficies refratoras
  PsiLC = [];         #Coeficiente de transicao, necessario para calcular o H de um usuario no room 2
  %==============================================================

  if usuario(1) > 0 # usuario esta no room 1
    cosseno_csi_AU = cosseno_csi(VLC, usuario, alpha, beta, d);
    cosseno_phi_AU = abs(usuario(3) - VLC(3)) / d;  #usando a formula que eu criei, podemos mudar para a que a giovanna encontrou para testar
    #Distância entre os espelhos e o usuário e entre espelhos e AP
    for i = 1:num_esp
      espelho = espelhos{i};
      dKU(i) = sqrt(sum((usuario - espelho).^2));
      dAK(i) = sqrt(sum((VLC - espelho).^2));
    endfor

    #calcula o cosseno_csi e o cosseno_phi para espelho
    for i = 1:num_esp
      espelho = espelhos{i};
      cossenos_phi_AK(i) = abs(espelho(3) - VLC(3)) / dAK(i);
      cossenos_csi_AK(i) = cosseno_csi(VLC, espelho, alpha, beta, dAK(i));
      cossenos_csi_KU(i) = cosseno_csi(espelho, usuario, alpha, beta, dKU(i));
      cossenos_phi_KU(i) = cosseno_phi(espelho, usuario, yaw, roll, dKU(i));
    endfor

    #Calcula HLos para o usuário
    HLos = Hlos(m, G_csi, T_csi, cosseno_phi_AU, cosseno_csi_AU, Apd, d);

    #HNlosk para cada espelho
    for i = 1:num_esp
      Hnlosk(i) = HNlosk(m, ro_RIS, Apd, dAK(i), dKU(i), d, Ak, G_csi, T_csi, cossenos_phi_AK(i), cossenos_csi_AK(i), cossenos_phi_KU(i), cossenos_csi_KU(i));
    endfor

    H = iota * HLos + sum(Hnlosk);
    Gamma = 0;
  else # usuario esta no room 2
    #Distância entre os refratores e o usuário e entre refratores e AP
    for i = 1:num_ref
      refrator = refratores{i};
      dNU(i) = sqrt(sum((usuario - refrator).^2));
      dAN(i) = sqrt(sum((VLC - refrator).^2));
    endfor

    for i = 1:num_ref
      refrator = refratores{i};
      cossenos_phi_AN(i) = abs(refrator(3) - VLC(3)) / dAN(i);
      cossenos_csi_AN(i) = cosseno_csi(VLC, refrator, alpha, beta, dAN(i));
      cossenos_csi_NU(i) = cosseno_csi(refrator, usuario, alpha, beta, dNU(i));
      senos_csi_NU(i) = sqrt(1 - cossenos_csi_NU(i)^2); #Identidade trigonométrica fundamental
      cossenos_phi_NU(i) = cosseno_phi(refrator, usuario, yaw, roll, dNU(i));
      senos_teta(i) = (senos_csi_NU(i) * Eta_a) / Eta_c; #Lei de snell
      cossenos_teta(i) = sqrt(1 - senos_teta(i)^2); #Identidade trigonométrica fundamental
    endfor

    #HNlosn para cada refrator
    for i = 1:num_ref
      Hnlosn(i) = HNlosn(m, Apd, dAN(i), dNU(i), d, An, G_csi, T_csi, cossenos_phi_AN(i), cossenos_csi_AN(i), cossenos_phi_NU(i), cossenos_csi_NU(i));
    endfor

    Eta = Eta_c/Eta_a;
    Eta1 = Eta_a/Eta_c;
    #calculo do PsiLc para cada refrator
    for i = 1:num_ref
      Tac(i) = (1-((((Eta^2 * cossenos_csi_NU(i) - sqrt(Eta^2 - senos_csi_NU(i)^2))/
              (Eta^2 * cossenos_csi_NU(i) + sqrt(Eta^2 - senos_csi_NU(i)^2)))^2)/2 +
              (((cossenos_csi_NU(i) - sqrt(Eta^2 - senos_csi_NU(i)^2))/
              (cossenos_csi_NU(i) + sqrt(Eta^2 - senos_csi_NU(i)^2)))^2)/2));
      Tca(i) = (1-((((cossenos_teta(i) - sqrt(Eta1^2 - senos_teta(i)^2))/
              (cossenos_teta(i) + sqrt(Eta1^2 - senos_teta(i)^2)))^2)/2 +
              (((Eta1^2 * cossenos_teta(i) - sqrt(Eta1^2 - senos_teta(i)^2))/
              (Eta1^2 * cossenos_teta(i) + sqrt(Eta1^2 - senos_teta(i)^2)))^2)/2));
      PsiLC(i) = Tac(i) * Tca(i);
    endfor

    H = sum(Hnlosn .* PsiLC);

    #calculo do Gamma para usuarios no room2, sera usado depois no calculo do sum rate
    aux = (Eta_o*sqrt((Eta_e^2-Eta_o^2)*(Eta_e^2-Eta_c^2))) / (Eta_c*(Eta_e^2-Eta_o^2));
    Ve = Vth - log10(-tan(tan(aux)^(-1)/2 - pi/4));
    E = Ve/D;
    Gamma = ((2*pi*Eta_c^2)/(cossenos_csi_NU(floor(num_ref/2))*lambda)) * r_eff*E;
    # o artigo usa cosseno_csi_NU, mas nao fala qual refrator N eh, entao estou meio confuso
  endif
endfunction

#Cálculo do cosseno do ângulo de incidência
function Cosseno_csi = cosseno_csi(coord1, coord2, alpha, beta, d)
  Cosseno_csi = ((abs(coord1(1) - coord2(1))/d) * cosd(beta) * sind(alpha) +
                (abs(coord1(2) - coord2(2))/d) * sind(beta) * sind(alpha) +
                (abs(coord1(3) - coord2(3))/d) * cosd(alpha));
endfunction

#Cálculo do cosseno do angulo de radiacao
function Cosseno_phi = cosseno_phi(coord1, coord2, yaw, roll, d)
  Cosseno_phi = (((abs(coord1(1) - coord2(1))/d)) * sind(yaw) * cosd(roll) +
                (abs(coord1(2) - coord2(2))/d) * cosd(yaw) * cosd(roll) +
                (abs(coord1(3) - coord2(3))/d)  *sind(roll));
endfunction

#Ganho do Canal para o Caminho da Linha de Visão (HLos)
#Consideramos que o ângulo de incidência csi está dentro do intervalo correto
# 0 <= Csi <= Csi_fov
#Caso não esteja, Hlos = 0
function DevolveHlos = Hlos(m, G_csi, T_csi, cosseno_phi, cosseno_csi, Apd, d)
  DevolveHlos = ((((m + 1) * Apd)/((2 * pi) * (d^2))) * G_csi * T_csi * (cosseno_phi)^m * cosseno_csi);
endfunction

#Ganho do Canal para o Caminho SEM Linha de Visão (HNlosk)
#Com superfícies refletoras
#Consideramos que o ângulo de incidência csi está dentro do intervalo correto
#G_csi e T_csi já estão definidos
function DevolveHNlosk = HNlosk(m, ro_RIS, Apd, dAK, dKU, d, Ak, G_csi, T_csi, cosseno_phi_AK, cosseno_csi_AK, cosseno_phi_KU, cosseno_csi_KU)
  DevolveHNlosk = (ro_RIS * (((m + 1) * Apd)/(2*(pi^2) * dAK^2 * dKU^2)) * d * Ak * G_csi * T_csi *
                  cosseno_phi_AK^m * cosseno_csi_AK * cosseno_phi_KU * cosseno_csi_KU);
endfunction

#Ganho do Canal SEM linha de visão (HNlosn)
#Ganho de canal para o para o caminho NLos para a propagação do sinal do AP
#através da n-ésima superfície REFRATORA e o usuario u
function DevolveHNlosn = HNlosn(m, Apd, dAN, dNU, d, An, G_csi, T_csi, cosseno_phi_AN, cosseno_csi_AN, cosseno_phi_NU, cosseno_csi_NU)
  DevolveHNlosn = ((((m + 1) * Apd)/(2*(pi^2) * dAN^2 * dNU^2)) * d * An * G_csi * T_csi *
                  cosseno_phi_AN^m * cosseno_csi_AN * cosseno_phi_NU * cosseno_csi_NU);
endfunction
