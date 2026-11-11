var Alineamiento = (function () {
  var abierto = -20;
  var extension = -1;
  var bruto = "AA:8,AR:-2,AN:-4,AD:-4,AC:0,AQ:-2,AE:-2,AG:0,AH:-4,AI:-2,AL:-2,AK:-2,AM:-2,AF:-4,AP:-2,AS:2,AT:0,AW:-6,AY:-4,AV:0,AB:-4,AZ:-2,AX:0,A*:-8,RA:-2,RR:10,RN:0,RD:-4,RC:-6,RQ:2,RE:0,RG:-4,RH:0,RI:-6,RL:-4,RK:4,RM:-2,RF:-6,RP:-4,RS:-2,RT:-2,RW:-6,RY:-4,RV:-6,RB:-2,RZ:0,RX:-2,R*:-8,NA:-4,NR:0,NN:12,ND:2,NC:-6,NQ:0,NE:0,NG:0,NH:2,NI:-6,NL:-6,NK:0,NM:-4,NF:-6,NP:-4,NS:2,NT:0,NW:-8,NY:-4,NV:-6,NB:6,NZ:0,NX:-2,N*:-8,DA:-4,DR:-4,DN:2,DD:12,DC:-6,DQ:0,DE:4,DG:-2,DH:-2,DI:-6,DL:-8,DK:-2,DM:-6,DF:-6,DP:-2,DS:0,DT:-2,DW:-8,DY:-6,DV:-6,DB:8,DZ:2,DX:-2,D*:-8,CA:0,CR:-6,CN:-6,CD:-6,CC:18,CQ:-6,CE:-8,CG:-6,CH:-6,CI:-2,CL:-2,CK:-6,CM:-2,CF:-4,CP:-6,CS:-2,CT:-2,CW:-4,CY:-4,CV:-2,CB:-6,CZ:-6,CX:-4,C*:-8,QA:-2,QR:2,QN:0,QD:0,QC:-6,QQ:10,QE:4,QG:-4,QH:0,QI:-6,QL:-4,QK:2,QM:0,QF:-6,QP:-2,QS:0,QT:-2,QW:-4,QY:-2,QV:-4,QB:0,QZ:6,QX:-2,Q*:-8,EA:-2,ER:0,EN:0,ED:4,EC:-8,EQ:4,EE:10,EG:-4,EH:0,EI:-6,EL:-6,EK:2,EM:-4,EF:-6,EP:-2,ES:0,ET:-2,EW:-6,EY:-4,EV:-4,EB:2,EZ:8,EX:-2,E*:-8,GA:0,GR:-4,GN:0,GD:-2,GC:-6,GQ:-4,GE:-4,GG:12,GH:-4,GI:-8,GL:-8,GK:-4,GM:-6,GF:-6,GP:-4,GS:0,GT:-4,GW:-4,GY:-6,GV:-6,GB:-2,GZ:-4,GX:-2,G*:-8,HA:-4,HR:0,HN:2,HD:-2,HC:-6,HQ:0,HE:0,HG:-4,HH:16,HI:-6,HL:-6,HK:-2,HM:-4,HF:-2,HP:-4,HS:-2,HT:-4,HW:-4,HY:4,HV:-6,HB:0,HZ:0,HX:-2,H*:-8,IA:-2,IR:-6,IN:-6,ID:-6,IC:-2,IQ:-6,IE:-6,IG:-8,IH:-6,II:8,IL:4,IK:-6,IM:2,IF:0,IP:-6,IS:-4,IT:-2,IW:-6,IY:-2,IV:6,IB:-6,IZ:-6,IX:-2,I*:-8,LA:-2,LR:-4,LN:-6,LD:-8,LC:-2,LQ:-4,LE:-6,LG:-8,LH:-6,LI:4,LL:8,LK:-4,LM:4,LF:0,LP:-6,LS:-4,LT:-2,LW:-4,LY:-2,LV:2,LB:-8,LZ:-6,LX:-2,L*:-8,KA:-2,KR:4,KN:0,KD:-2,KC:-6,KQ:2,KE:2,KG:-4,KH:-2,KI:-6,KL:-4,KK:10,KM:-2,KF:-6,KP:-2,KS:0,KT:-2,KW:-6,KY:-4,KV:-4,KB:0,KZ:2,KX:-2,K*:-8,MA:-2,MR:-2,MN:-4,MD:-6,MC:-2,MQ:0,ME:-4,MG:-6,MH:-4,MI:2,ML:4,MK:-2,MM:10,MF:0,MP:-4,MS:-2,MT:-2,MW:-2,MY:-2,MV:2,MB:-6,MZ:-2,MX:-2,M*:-8,FA:-4,FR:-6,FN:-6,FD:-6,FC:-4,FQ:-6,FE:-6,FG:-6,FH:-2,FI:0,FL:0,FK:-6,FM:0,FF:12,FP:-8,FS:-4,FT:-4,FW:2,FY:6,FV:-2,FB:-6,FZ:-6,FX:-2,F*:-8,PA:-2,PR:-4,PN:-4,PD:-2,PC:-6,PQ:-2,PE:-2,PG:-4,PH:-4,PI:-6,PL:-6,PK:-2,PM:-4,PF:-8,PP:14,PS:-2,PT:-2,PW:-8,PY:-6,PV:-4,PB:-4,PZ:-2,PX:-4,P*:-8,SA:2,SR:-2,SN:2,SD:0,SC:-2,SQ:0,SE:0,SG:0,SH:-2,SI:-4,SL:-4,SK:0,SM:-2,SF:-4,SP:-2,SS:8,ST:2,SW:-6,SY:-4,SV:-4,SB:0,SZ:0,SX:0,S*:-8,TA:0,TR:-2,TN:0,TD:-2,TC:-2,TQ:-2,TE:-2,TG:-4,TH:-4,TI:-2,TL:-2,TK:-2,TM:-2,TF:-4,TP:-2,TS:2,TT:10,TW:-4,TY:-4,TV:0,TB:-2,TZ:-2,TX:0,T*:-8,WA:-6,WR:-6,WN:-8,WD:-8,WC:-4,WQ:-4,WE:-6,WG:-4,WH:-4,WI:-6,WL:-4,WK:-6,WM:-2,WF:2,WP:-8,WS:-6,WT:-4,WW:22,WY:4,WV:-6,WB:-8,WZ:-6,WX:-4,W*:-8,YA:-4,YR:-4,YN:-4,YD:-6,YC:-4,YQ:-2,YE:-4,YG:-6,YH:4,YI:-2,YL:-2,YK:-4,YM:-2,YF:6,YP:-6,YS:-4,YT:-4,YW:4,YY:14,YV:-2,YB:-6,YZ:-4,YX:-2,Y*:-8,VA:0,VR:-6,VN:-6,VD:-6,VC:-2,VQ:-4,VE:-4,VG:-6,VH:-6,VI:6,VL:2,VK:-4,VM:2,VF:-2,VP:-4,VS:-4,VT:0,VW:-6,VY:-2,VV:8,VB:-6,VZ:-4,VX:-2,V*:-8,BA:-4,BR:-2,BN:6,BD:8,BC:-6,BQ:0,BE:2,BG:-2,BH:0,BI:-6,BL:-8,BK:0,BM:-6,BF:-6,BP:-4,BS:0,BT:-2,BW:-8,BY:-6,BV:-6,BB:8,BZ:2,BX:-2,B*:-8,ZA:-2,ZR:0,ZN:0,ZD:2,ZC:-6,ZQ:6,ZE:8,ZG:-4,ZH:0,ZI:-6,ZL:-6,ZK:2,ZM:-2,ZF:-6,ZP:-2,ZS:0,ZT:-2,ZW:-6,ZY:-4,ZV:-4,ZB:2,ZZ:8,ZX:-2,Z*:-8,XA:0,XR:-2,XN:-2,XD:-2,XC:-4,XQ:-2,XE:-2,XG:-2,XH:-2,XI:-2,XL:-2,XK:-2,XM:-2,XF:-2,XP:-4,XS:0,XT:0,XW:-4,XY:-2,XV:-2,XB:-2,XZ:-2,XX:-2,X*:-8,*A:-8,*R:-8,*N:-8,*D:-8,*C:-8,*Q:-8,*E:-8,*G:-8,*H:-8,*I:-8,*L:-8,*K:-8,*M:-8,*F:-8,*P:-8,*S:-8,*T:-8,*W:-8,*Y:-8,*V:-8,*B:-8,*Z:-8,*X:-8,**:2";
  var matriz = {};
  bruto.split(",").forEach(function (par) {
    matriz[par.slice(0, 2)] = parseInt(par.slice(3), 10);
  });
  var negativo = -1000000000;
  function puntuar(a, b) {
    var n = a.length;
    var m = b.length;
    var M = [];
    var X = [];
    var Y = [];
    for (var i = 0; i <= n; i++) {
      M.push(new Array(m + 1).fill(negativo));
      X.push(new Array(m + 1).fill(negativo));
      Y.push(new Array(m + 1).fill(negativo));
    }
    M[0][0] = 0;
    for (i = 1; i <= n; i++) X[i][0] = abierto + (i - 1) * extension;
    for (var j = 1; j <= m; j++) Y[0][j] = abierto + (j - 1) * extension;
    for (i = 1; i <= n; i++) {
      for (j = 1; j <= m; j++) {
        var s = matriz[a[i - 1] + b[j - 1]];
        if (s === undefined) s = negativo;
        X[i][j] = Math.max(M[i - 1][j] + abierto, X[i - 1][j] + extension, Y[i - 1][j] + abierto);
        Y[i][j] = Math.max(M[i][j - 1] + abierto, Y[i][j - 1] + extension, X[i][j - 1] + abierto);
        M[i][j] = s + Math.max(M[i - 1][j - 1], X[i - 1][j - 1], Y[i - 1][j - 1]);
      }
    }
    return { M: M, X: X, Y: Y };
  }
  function alinear(a, b) {
    var tablas = puntuar(a, b);
    var M = tablas.M;
    var X = tablas.X;
    var Y = tablas.Y;
    var i = a.length;
    var j = b.length;
    var estados = { M: M, X: X, Y: Y };
    var estado = "M";
    if (X[i][j] > estados[estado][i][j]) estado = "X";
    if (Y[i][j] > estados[estado][i][j]) estado = "Y";
    var sa = [];
    var sb = [];
    while (i || j) {
      if (estado === "M") {
        sa.push(a[i - 1]);
        sb.push(b[j - 1]);
        i -= 1;
        j -= 1;
        if (i || j) {
          estado = "M";
          if (X[i][j] > estados[estado][i][j]) estado = "X";
          if (Y[i][j] > estados[estado][i][j]) estado = "Y";
        }
      } else if (estado === "X") {
        sa.push(a[i - 1]);
        sb.push("-");
        var siguiente = "M";
        if (M[i - 1][j] + abierto !== X[i][j]) {
          if (X[i - 1][j] + extension === X[i][j]) siguiente = "X";
          else siguiente = "Y";
        }
        estado = siguiente;
        i -= 1;
      } else {
        sa.push("-");
        sb.push(b[j - 1]);
        siguiente = "M";
        if (M[i][j - 1] + abierto !== Y[i][j]) {
          if (Y[i][j - 1] + extension === Y[i][j]) siguiente = "Y";
          else siguiente = "X";
        }
        estado = siguiente;
        j -= 1;
      }
    }
    sa.reverse();
    sb.reverse();
    var identidades = 0;
    var similares = 0;
    var huecos = 0;
    for (var k = 0; k < sa.length; k++) {
      if (sa[k] === "-" || sb[k] === "-") huecos += 1;
      else if (sa[k] === sb[k]) identidades += 1;
      if (sa[k] !== "-" && sb[k] !== "-" && matriz[sa[k] + sb[k]] > 0) similares += 1;
    }
    var mejor = Math.max(M[a.length][b.length], X[a.length][b.length], Y[a.length][b.length]);
    return {
      puntuacion: mejor / 2,
      identidades: identidades,
      similares: similares,
      huecos: huecos,
      columnas: sa.length,
      arriba: sa.join(""),
      abajo: sb.join("")
    };
  }
  return { alinear: alinear };
})();
