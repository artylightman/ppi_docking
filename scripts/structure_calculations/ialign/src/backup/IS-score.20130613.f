******************************************************************************
*     This program compares two protein-protein interfaces and identifies the 
*     best superimposition that has the maximum iTM-score or IS-score. It is  
*     based on TM-score originally written by Yang Zhang. You may modify and
*     distribute the source code for academic research purpose only.
*
*     Please address your comments and question to: Mu Gao <mu.gao@gatech.edu>
*
*     References: 
*     Mu Gao, Jeffrey Skolnick, Proteins, 2011 79:1623-34.
*     Mu Gao, Jeffrey Skolnick, Bioinformatics, 2010 26:2259-65.
*     Yang Zhang, Jeffrey Skolnick, Proteins, 2004 57:702-10.
*
******************* Updating history ************************************
*
*     2010/07/21: initial distribution
*     2010/12/10: various bug fix
*************************************************************************
      
      program ISscore
      PARAMETER(nmax=3000)
      PARAMETER(cmax=100)   !Mu: maximum number of contact per residue
      PARAMETER(maxk=20)    !Mu: maximum number of chains
      
      common/stru/xt(nmax),yt(nmax),zt(nmax),xb(nmax),yb(nmax),zb(nmax)
      common/nres/nseqA,nseqB
      common/para/d,d0,d0_fix
      common/align/n_ali,iA(nmax),iB(nmax)
      common/nscore/i_ali(nmax),n_cut ![1,n_ali],align residues for the score
      dimension k_ali(nmax),k_ali0(nmax)

      character*500 fnam,pdb(100),outname,file1,file2,matchfile
      character*3 aa(-1:20),seqA(nmax),seqB(nmax)
      character*5 nresA(nmax),nresB(nmax)
      character*100 s,du
      character seq1A(nmax),seq1B(nmax),ali(nmax),chA(nmax),chB(nmax)
      character sequenceA(nmax),sequenceB(nmax),sequenceM(nmax)
      character*2   measure   ! 'TM' -> TM-score, 'IS' -> IS-score
      character*80  matchlst(nmax)    !list of matched residues
      character   pdbchnmA(maxk), pdbchnmB(maxk)

      character*5   matresid1(nmax),matresid2(nmax)
      character     matresch1(nmax),matresch2(nmax)

      dimension L_ini(100),iq(nmax)
      common/scores/score,score_maxsub,score_fix,score10
      common/GDT/n_GDT05,n_GDT1,n_GDT2,n_GDT4,n_GDT8
      double precision score,score_max,score_fix,score_fix_max
      double precision score_maxsub,score10
      dimension xa(nmax),ya(nmax),za(nmax)

ccc   Contact stuff added by Mu
      common/contcf/d_col
      common/contol/numcol                !total number of overlapped contacts
      common/id2chain/id_chain1(nmax),id_chain2(nmax)  !original index to chain ID
      integer cont1(cmax,nmax),cont2(cmax,nmax)     !contact list
      common/contlst/cont1,cont2
      common/contnum/ncont1(nmax),ncont2(nmax)      !number of contacts
      character*100 contfile1,contfile2    !Mu: contact list files
      integer numcol_max

      common/chains/ichtermA(0:maxk),ichtermB(0:maxk) !indexes of C-termini of chains

ccc   RMSD:
      double precision r_1(3,nmax),r_2(3,nmax),r_3(3,nmax),w(nmax)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /nmax*1.0/
ccc   

      data aa/ 'BCK','GLY','ALA','SER','CYS',
     &     'VAL','THR','ILE','PRO','MET',
     &     'ASP','ASN','LEU','LYS','GLU',
     &     'GLN','ARG','HIS','PHE','TYR',
     &     'TRP','CYX'/
      character*1 slc(-1:20)
      data slc/'X','G','A','S','C',
     &     'V','T','I','P','M',
     &     'D','N','L','K','E',
     &     'Q','R','H','F','Y',
     &     'W','C'/

ccc   Options:
      integer  score_flag  ! 0 - TM-score, 1 - IS-score
      integer  eqpos_flag  ! 0 - identical residue id are equivilent, 1 - residues with the same index are equivilent
      integer  verbo_flag  ! 2 - print detailed match list, default 1
      integer  match_flag  ! 1 - read a list of matched residues, default 0

      call getarg(1,fnam)
      if(fnam.eq.' '.or.fnam.eq.'?'.or.fnam.eq.'-h')then
         call printHelp()
         goto 9999
      endif
      
******* options ----------->
      m_out=-1
      m_fix=-1
      score_flag=1        !default: 1 - IS-score
      eqpos_flag=0        !default: 0 - residues with same identity are equivilent
      verbo_flag=1
      match_flag=0
      measure='IS'
      narg=iargc()
      i=0
      j=0
 115  continue
      i=i+1
      call getarg(i,fnam)
      if(fnam.eq.'-d')then
         m_fix=1
         i=i+1
         call getarg(i,fnam)
         read(fnam,*)d0_fix
      elseif(fnam.eq.'-t')then ! use regular TM-score measure
         score_flag=0
         measure='TM'
      elseif(fnam.eq.'-eqpos')then 
         eqpos_flag=1
      elseif(fnam.eq.'-v')then ! print detailed match list
         i=i+1
         call getarg(i,fnam)
         read(fnam,*)verbo_flag
         if(verbo_flag.lt.0.or.verbo_flag.gt.2)then
            write(*,*)'Eorror: invalid verbose option (valid: 0,1,2)'
            goto 9999
         endif
      elseif(fnam.eq.'-m')then
         i=i+1
         match_flag=1
         call getarg(i,matchfile)
      else
         j=j+1
         pdb(j)=fnam
      endif
      if(i.lt.narg)goto 115

ccccccccc read data from first CA file:
      open(unit=10,file=pdb(1),status='old')
      i=0
      ich=0
      ichtermA(ich)=0
 101  read(10,104,end=102) s
c      if(s(1:3).eq.'TER') goto 102   ! Mu: consider multi-chains
      if(s(1:4).eq.'ATOM')then
         if(s(13:16).eq.'CA  '.or.s(13:16).eq.' CA '.or.s(13:16).
     &        eq.'  CA')then
         if(s(17:17).eq.' '.or.s(17:17).eq.'A')then
            i=i+1
            read(s,103)du,seqA(i),du,chA(i),nresA(i),du,xa(i),ya(i),za(i)

            if(i.gt.1)then
               if(chA(i).ne.chA(i-1))then
                  ich=ich+1
                  ichtermA(ich)=i-1
                  pdbchnmA(ich)=chA(i-1)
               endif
            endif

            do j=-1,20
               if(seqA(i).eq.aa(j))then
                  seq1A(i)=slc(j)
                  goto 21
               endif
            enddo
            seq1A(i)=slc(-1)
 21         continue
         endif
         endif
      endif
      goto 101
 102  continue
 103  format(A17,A3,A1,A1,A5,A3,3F8.3)
 104  format(A100)
      close(10)
      nseqA=i
      ichtermA(ich+1)=i
      pdbchnmA(ich+1)=chA(i)
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      
ccccccccc read data from first CA file:
      open(unit=10,file=pdb(2),status='old')
      i=0
      ich=0
      ichtermB(ich)=0
 201  read(10,204,end=202) s
c      if(s(1:3).eq.'TER') goto 202   ! Mu: consider multi-chains
      if(s(1:4).eq.'ATOM')then
         if(s(13:16).eq.'CA  '.or.s(13:16).eq.' CA '.or.s(13:16).
     &        eq.'  CA')then
         if(s(17:17).eq.' '.or.s(17:17).eq.'A')then
            i=i+1
            read(s,203)du,seqB(i),du,chB(i),nresB(i),du,xb(i),yb(i),zb(i)

            if(i.gt.1)then
               if(chB(i).ne.chB(i-1))then
                  ich=ich+1
                  ichtermB(ich)=i-1
                  pdbchnmB(ich)=chB(i-1)
               endif
            endif

            do j=-1,20
               if(seqB(i).eq.aa(j))then
                  seq1B(i)=slc(j)
                  goto 22
               endif
            enddo
            seq1B(i)=slc(-1)
 22         continue
         endif
         endif
      endif
      goto 201
 202  continue
 203  format(A17,A3,A1,A1,A5,A3,3F8.3)
 204  format(A100)
      close(10)
      nseqB=i
      ichtermB(ich+1)=i
      pdbchnmB(ich+1)=chB(i)
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


      contfile1=pdb(3)
      contfile2=pdb(4)
      call readContact(contfile1,contfile2)   !Mu: read contact lists pre-calculated

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      num_cont1=0  !number of interfacial contacts of protein 1
      num_cont2=0  !number of interfacial contacts of protein 2
      do i=1,nseqA
         if(i.le.ichtermA(1))then
            id_chain1(i)=1
            num_cont1=num_cont1 + ncont1(i)  ! only need count once
         else
            id_chain1(i)=2
         endif
      enddo
      do j=1,nseqB
         if(j.le.ichtermB(1))then
            id_chain2(j)=1
            num_cont2=num_cont2 + ncont2(j)
         else
            id_chain2(j)=2
         endif
      enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




******************************************************************
*     pickup the aligned residues:
******************************************************************
      if(match_flag.eq.1)then
         call readMatchFile(matchfile,matresch1,matresid1,
     &        matresch2,matresid2,nmatch)
      endif

      k=0
      if(eqpos_flag.eq.0)then
         if(match_flag.eq.1)then  !use the match list provided
         do kk=1,nmatch
            do i=1,nseqA
               if(nresA(i).eq.matresid1(kk).and.chA(i).eq.matresch1(kk))then  
                  iA(kk)=i
                  goto 205
               endif
            enddo
 205        continue

            do j=1,nseqB
               if(nresB(j).eq.matresid2(kk).and.chB(j).eq.matresch2(kk))then  
                  iB(kk)=j
                  goto 206
               endif
            enddo
 206        continue
         enddo
         k=nmatch

         else  ! match list not provided, residues with the same resid id and chain id
         do i=1,nseqA
            do j=1,nseqB
               if(nresA(i).eq.nresB(j).and.chA(i).eq.chB(j))then  
                  k=k+1
                  iA(k)=i
                  iB(k)=j
                  goto 207
               endif
            enddo
 207        continue
         enddo
         endif
      else
c     Mu: residues at the same position (index starting from 0) are equivilent
         if( nseqA .le. nseqB ) then
            do i=1,nseqA
               k=k+1
               iA(k)=i
               iB(k)=i
            enddo
         else
            do i=1,nseqB
               k=k+1
               iA(k)=i
               iB(k)=i
            enddo
         endif
      endif

      n_ali=k                   !number of aligned residues
      if(n_ali.lt.1)then
        write(*,*)'There is no common residue in the input structures'
        goto 9999
      endif


************/////
*     parameters:
*****************
***   d0------------->
      if(nseqB.gt.15)then
         d0=1.24*(nseqB-15)**(1.0/3.0)-1.8
      else
         d0=0.5
      endif
      if(d0.lt.0.5)d0=0.5
      if(m_fix.eq.1)d0=d0_fix
***   d0_search ----->
      d0_search=d0
      if(d0_search.gt.8)d0_search=8
      if(d0_search.lt.4.5)d0_search=4.5

      d_col=d0_search  !distance cutoff for identifying overlaped contact

***   iterative parameters ----->
      n_it=20                   !maximum number of iterations
      d_output=5                !for output alignment
      if(m_fix.eq.1)d_output=d0_fix
      n_init_max=6              !maximum number of L_init
      n_init=0
      L_ini_min=4
      if(n_ali.lt.4)L_ini_min=n_ali
      do i=1,n_init_max-1
         n_init=n_init+1
         L_ini(n_init)=n_ali/2**(n_init-1)
         if(L_ini(n_init).le.L_ini_min)then
            L_ini(n_init)=L_ini_min
            goto 402
         endif
      enddo
      n_init=n_init+1
      L_ini(n_init)=L_ini_min
 402  continue
      

******************************************************************
*     find all interfacial contacts
******************************************************************



******************************************************************
*     find the maximum score starting from local structures superposition
******************************************************************
      score_max=-1              !TM/IS-score
      score_maxsub_max=-1       !MaxSub-score
      score10_max=-1            !TM-score10
      numcol_max=-1             !contact overlap
      n_GDT05_max=-1            !number of residues<0.5
      n_GDT1_max=-1             !number of residues<1
      n_GDT2_max=-1             !number of residues<2
      n_GDT4_max=-1             !number of residues<4
      n_GDT8_max=-1             !number of residues<8
      do 333 i_init=1,n_init
        L_init=L_ini(i_init)
        iL_max=n_ali-L_init+1
        do 300 iL=1,iL_max      !on aligned residues, [1,nseqA]
           LL=0
           ka=0
           do i=1,L_init
              k=iL+i-1          ![1,n_ali] common aligned
              r_1(1,i)=xa(iA(k))
              r_1(2,i)=ya(iA(k))
              r_1(3,i)=za(iA(k))
              r_2(1,i)=xb(iB(k))
              r_2(2,i)=yb(iB(k))
              r_2(3,i)=zb(iB(k))
              ka=ka+1
              k_ali(ka)=k
              LL=LL+1
           enddo
           call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2
           if(i_init.eq.1)then  !global superposition
              armsd=dsqrt(rms/LL)
              rmsd_ali=armsd
           endif
           do j=1,nseqA
              xt(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
              yt(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
              zt(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
           enddo
           d=d0_search-1
           call score_fun (score_flag)      !init, get scores, n_cut+i_ali(i) for iteration
           if(score_max.lt.score)then
              score_max=score
              numcol_max=numcol
              ka0=ka
              do i=1,ka0
                 k_ali0(i)=k_ali(i)
              enddo
           endif
           if(score10_max.lt.score10)score10_max=score10
           if(score_maxsub_max.lt.score_maxsub)score_maxsub_max=
     &          score_maxsub
           if(n_GDT05_max.lt.n_GDT05)n_GDT05_max=n_GDT05
           if(n_GDT1_max.lt.n_GDT1)n_GDT1_max=n_GDT1
           if(n_GDT2_max.lt.n_GDT2)n_GDT2_max=n_GDT2
           if(n_GDT4_max.lt.n_GDT4)n_GDT4_max=n_GDT4
           if(n_GDT8_max.lt.n_GDT8)n_GDT8_max=n_GDT8

***   iteration for extending ---------------------------------->
           d=d0_search+1
           do 301 it=1,n_it
              LL=0
              ka=0
              do i=1,n_cut
                 m=i_ali(i)     ![1,n_ali]
                 r_1(1,i)=xa(iA(m))
                 r_1(2,i)=ya(iA(m))
                 r_1(3,i)=za(iA(m))
                 r_2(1,i)=xb(iB(m))
                 r_2(2,i)=yb(iB(m))
                 r_2(3,i)=zb(iB(m))
                 ka=ka+1
                 k_ali(ka)=m
                 LL=LL+1
              enddo
              call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2
              do j=1,nseqA
                 xt(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
                 yt(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
                 zt(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
              enddo
              call score_fun(score_flag)    !get scores, n_cut+i_ali(i) for iteration
              if(score_max.lt.score)then
                 score_max=score
                 numcol_max=numcol
                 ka0=ka
                 do i=1,ka
                    k_ali0(i)=k_ali(i)
                 enddo
              endif
              if(score10_max.lt.score10)score10_max=score10
              if(score_maxsub_max.lt.score_maxsub)score_maxsub_max
     &             =score_maxsub
              if(n_GDT05_max.lt.n_GDT05)n_GDT05_max=n_GDT05
              if(n_GDT1_max.lt.n_GDT1)n_GDT1_max=n_GDT1
              if(n_GDT2_max.lt.n_GDT2)n_GDT2_max=n_GDT2
              if(n_GDT4_max.lt.n_GDT4)n_GDT4_max=n_GDT4
              if(n_GDT8_max.lt.n_GDT8)n_GDT8_max=n_GDT8

 303          format(i5,i5,i5,f17.14,f17.14,i5,f7.3)
              if(it.eq.n_it)goto 302
              if(n_cut.eq.ka)then
                 neq=0
                 do i=1,n_cut
                    if(i_ali(i).eq.k_ali(i))neq=neq+1
                 enddo
                 if(n_cut.eq.neq)goto 302
              endif
 301       continue             !for iteration
 302       continue
 300    continue                !for shift
 333  continue                  !for initial length, L_ali/M

      if(score_flag.eq.1)then
         !rescale IS-score
         !f0=0.18 - 0.35*(nseqB**(-0.3))
	 f0=0.14-0.2*(nseqB**(-0.3))
         score_max=(score_max+f0)/(1+f0) 
         call calcPvalue(score_max,nseqB,pvalue,zs)
      endif


******************************************************************
*     Output
******************************************************************
***   output TM-scale ---------------------------->
      write(*,*)
      write(*,*)'*****************************************************',
     &     '************************'
      write(*,*)'                      Interface Similarity Score     '
      write(*,*)'  A metric for assessing the quality of protein',
     &     ' protein interface models.'
      write(*,*)
      write(*,*)'  Reference: Mu Gao and Jeffrey Skolnick, ',
     &     'Bioinformatics 2010, in press'
      write(*,*)'  Please address your omments to <mu.gao@gatech.edu>.'
      write(*,*)'*****************************************************',
     &     '************************'
      write(*,*)

      call getBasename(pdb(1),file1)
      call getBasename(pdb(2),file2)

      write(*,501)1,file1,pdbchnmA(1),pdbchnmA(2),nseqA,num_cont1
      write(*,501)2,file2,pdbchnmB(1),pdbchnmB(2),nseqB,num_cont2
 501  format('Structure ',I1,': ',A25,' Chains ',A1,1X,A1,','
     &     I4, ' AAs,',I4,' Contacts')
      write(*,503)n_ali
 503  format('Number of matched residues = ',I4)
      write(*,513)rmsd_ali
 513  format('RMSD   of matched residues = ',F8.3)
      write(*,514)numcol_max
 514  format('Number of matched contacts = ',I4)
      write(*,*)
      if(score_flag.eq.1)then
        write(*,509)measure,score_max,pvalue,zs
      else
        write(*,510)measure,score_max
      endif
 509  format(a2,'-score    = ',f6.4, 
     &     '  P-value =',E12.4E3, '  Z-score =',f7.3)
 510  format(a2,'-score    = ',f6.4)
c      write(*,504)score_max,d0,score10_max
c 504  format('TM-score    = ',f6.4,'  (d0=',f5.2,',',' TM10= ',f6.4,')')
      write(*,505)score_maxsub_max
 505  format('MaxSub-score= ',f6.4,'  (d0= 3.50)')
      score_GDT=(n_GDT1_max+n_GDT2_max+n_GDT4_max+n_GDT8_max)
     &     /float(4*nseqB)
      write(*,506)score_GDT,n_GDT1_max/float(nseqB),n_GDT2_max
     &     /float(nseqB),n_GDT4_max/float(nseqB),n_GDT8_max/float(nseqB)
 506  format('GDT-TS-score= ',f6.4,' %(d<1)=',f6.4,' %(d<2)=',f6.4,
     $     ' %(d<4)=',f6.4,' %(d<8)=',f6.4)
      score_GDT_HA=(n_GDT05_max+n_GDT1_max+n_GDT2_max+n_GDT4_max)
     &     /float(4*nseqB)
      write(*,507)score_GDT_HA,n_GDT05_max/float(nseqB),n_GDT1_max
     &     /float(nseqB),n_GDT2_max/float(nseqB),n_GDT4_max/float(nseqB)
 507  format('GDT-HA-score= ',f6.4,' %(d<0.5)=',f6.4,' %(d<1)=',f6.4,
     $     ' %(d<2)=',f6.4,' %(d<4)=',f6.4)
      write(*,*)
      
***   recall and output the superposition of maxiumum TM-score:
      LL=0
      do i=1,ka0
         m=k_ali0(i)            !record of the best alignment
         r_1(1,i)=xa(iA(m))
         r_1(2,i)=ya(iA(m))
         r_1(3,i)=za(iA(m))
         r_2(1,i)=xb(iB(m))
         r_2(2,i)=yb(iB(m))
         r_2(3,i)=zb(iB(m))
         LL=LL+1
      enddo
      call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2
      do j=1,nseqA
         xt(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
         yt(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
         zt(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
      enddo

********* extract rotation matrix ------------>
      write(*,*)'-------- rotation matrix to rotate Chain-1 to ',
     &     'Chain-2 ------'
      write(*,*)'i          t(i)         u(i,1)         u(i,2) ',
     &     '        u(i,3)'
      do i=1,3
         write(*,304)i,t(i),u(i,1),u(i,2),u(i,3)
      enddo
c      do j=1,nseqA
c         xt(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
c         yt(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
c         zt(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
c         write(*,*)j,xt(j),yt(j),zt(j)
c      enddo
      write(*,*)
 304  format(I2,f18.10,f15.10,f15.10,f15.10)

********* rmsd in superposed regions --------------->
      d=d_output                !for output
      call score_fun(score_flag)          !give i_ali(i), score_max=score now
      LL=0
      do i=1,n_cut
         m=i_ali(i)             ![1,nseqA]
         r_1(1,i)=xa(iA(m))
         r_1(2,i)=ya(iA(m))
         r_1(3,i)=za(iA(m))
         r_2(1,i)=xb(iB(m))
         r_2(2,i)=yb(iB(m))
         r_2(3,i)=zb(iB(m))
         LL=LL+1
      enddo
      call u3b(w,r_1,r_2,LL,0,rms,u,t,ier)
      armsd=dsqrt(rms/LL)
      rmsd=armsd
      



*******************************************************************
***   output aligned sequences

      if(verbo_flag.eq.2) then
         call get_matchlst(nseqA,nseqB,iA,iB,xt,xb,yt,yb,zt,zb,
     &        nresA,nresB,seqA,seqB,ncont1,ncont2,chA,chB,n_ali,
     &        d_output,matchlst)

         write(*,*)'-----   Aligned Interface Residues  -----'
         write(*,*)'----- Structure 1       Structure 2 -----'
         write(*,*)'Index Ch1 Resid1 AA1    Ch2 Resid2 AA2   ',
     &     ' Distance NAC NC1 NC2 Note'

         do i = 1, nseqB
            write(*,'(A70)') matchlst(i)
         enddo

         write(*,*)
         write(*,'(A,F3.1,A)')'":" AA pairs within ',d_output,
     &        ' Ansgtrom in Ca-Ca distance'
         write(*,'(A)')'"*" Identical AA pairs, "-" Unmatched AA.'

      else if(verbo_flag.eq.1) then
         write(*,603) d_output
 603     format('Interface Alignment (":" AA pairs within ',F3.1,
     &        ' Angstrom in Ca-Ca distance)')
         write(*,608)pdbchnmA(1),pdbchnmA(2),pdbchnmB(1),pdbchnmB(2)
 608     format('Struct 1 Chains ',A1,A1,1X,' vs. ',
     &        'Struct 2 Chains ',A1,A1,
     &        ' (AAs in upper/lower cases)',/)

         call get_alnseq(iA,iB,n_ali,seq1A,seq1B,2,2,
     &        ichtermA,ichtermB,xt,yt,zt,xb,yb,zb,
     &        d_output,sequenceA,sequenceB,sequenceM,lseq)

         write(*,601)(sequenceA(i),i=1,lseq)
         write(*,601)(sequenceM(i),i=1,lseq)
         write(*,601)(sequenceB(i),i=1,lseq)
 601     format(2000A1)
      endif

      write(*,*)

      if(score_flag.eq.1)then
         write(*,109)nseqB,d0,f0
      else
         write(*,110)nseqB,d0
      endif

 109  format('Scoring parameters: normalization length =',I4,
     &     ',  d0 =',F6.3, ',  s0 =',F6.3)
 110  format('Scoring parameters: normalization length =',I4,
     &     ',  d0 =',F6.3)
      write(*,600)measure,n_cut,rmsd
 600  format('Superposed residues for the ',a2,'-score: length=',
     $     i3,', RMSD=',f6.2)
      write(*,*)

c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 9999 END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     1, collect those residues with dis<d;
c     2, calculate score_GDT, score_maxsub, score_TM
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine score_fun(score_flag)
      PARAMETER(nmax=3000)

      common/stru/xt(nmax),yt(nmax),zt(nmax),xb(nmax),yb(nmax),zb(nmax)
      common/nres/nseqA,nseqB
      common/para/d,d0,d0_fix
      common/align/n_ali,iA(nmax),iB(nmax)
      common/nscore/i_ali(nmax),n_cut ![1,n_ali],align residues for the score
      common/scores/score,score_maxsub,score_fix,score10
      common/GDT/n_GDT05,n_GDT1,n_GDT2,n_GDT4,n_GDT8
      double precision score,score_max,score_fix,score_fix_max
      double precision score_maxsub,score10

      dimension dis(nmax)

ccc   Contact stuff added by Mu
      common/contol/numcol                !total number of overlapped contacts
      common/contcf/d_col
      real fcol                           !contact overlap factor
      integer imap1(nmax),imap2(nmax)     !aligned position to original position
      integer is_aliA(nmax),is_aliB(nmax) !original position to aligned position
      integer col

      integer score_flag

      call init_ali(is_aliA,is_aliB,nseqA,nseqB)

      d_tmp=d
 21   n_cut=0                   !number of residue-pairs dis<d, for iteration
      n_GDT05=0                 !for GDT-score, # of dis<0.5
      n_GDT1=0                  !for GDT-score, # of dis<1
      n_GDT2=0                  !for GDT-score, # of dis<2
      n_GDT4=0                  !for GDT-score, # of dis<4
      n_GDT8=0                  !for GDT-score, # of dis<8
      score_maxsub_sum=0        !Maxsub-score
      score_sum=0               !TMscore
      score_sum10=0             !TMscore10
      do k=1,n_ali
         i=iA(k)                ![1,nseqA] reoder number of structureA
         j=iB(k)                ![1,nseqB]
         dis(k)=sqrt((xt(i)-xb(j))**2+(yt(i)-yb(j))**2+(zt(i)-zb(j))**2)
***   for iteration:
         if(dis(k).lt.d_tmp)then
            n_cut=n_cut+1
            i_ali(n_cut)=k      ![1,n_ali], mark the residue-pairs in dis<d
         endif
***   for GDT-score:
         if(dis(k).le.8)then
            n_GDT8=n_GDT8+1
            if(dis(k).le.4)then
               n_GDT4=n_GDT4+1
               if(dis(k).le.2)then
                  n_GDT2=n_GDT2+1
                  if(dis(k).le.1)then
                     n_GDT1=n_GDT1+1
                     if(dis(k).le.0.5)then
                        n_GDT05=n_GDT05+1
                     endif
                  endif
               endif
            endif
         endif
***   for MAXsub-score:
         if(dis(k).lt.3.5)then
            score_maxsub_sum=score_maxsub_sum+1/(1+(dis(k)/3.5)**2)
         endif
***   for TM-score:
         score_sum=score_sum+1/(1+(dis(k)/d0)**2)
***   for TM-score10:
         if(dis(k).lt.10)then
            score_sum10=score_sum10+1/(1+(dis(k)/d0)**2)
         endif
***   for IS-score:
         !if(dis(k).lt.d_col)then
            is_aliA(i)=k
            is_aliB(j)=k
         !endif
      enddo
      if(n_cut.lt.3.and.n_ali.gt.3)then
         d_tmp=d_tmp+.5
         goto 21
      endif

      score_is=0
      numcol=0
      do k=1,n_ali
         i=iA(k)
         j=iB(k)
***   calculate contact overlap
         call get_fcol(k,k,is_aliA,is_aliB,iA,iB,fcol,col)
***   IS-score
         score_is=score_is+fcol/(1+(dis(k)/d0)**2)
         numcol=numcol+col
      enddo

      score_maxsub=score_maxsub_sum/float(nseqB) !MAXsub-score
      if(score_flag.eq.0)then
         score=score_sum/real(nseqB) !TM-score
         score10=score_sum10/float(nseqB) !TM-score10
      else
         score=score_is/real(nseqB) !IS-score
      endif
      numcol=numcol/2
      
      return
      end

cccccccccccccccc Calculate sum of (r_d-r_m)^2 cccccccccccccccccccccccccc
c  w    - w(m) is weight for atom pair  c m           (given)
c  x    - x(i,m) are coordinates of atom c m in set x       (given)
c  y    - y(i,m) are coordinates of atom c m in set y       (given)
c  n    - n is number of atom pairs                         (given)
c  mode  - 0:calculate rms only                             (given)
c          1:calculate rms,u,t                              (takes longer)
c  rms   - sum of w*(ux+t-y)**2 over all atom pairs         (result)
c  u    - u(i,j) is   rotation  matrix for best superposition  (result)
c  t    - t(i)   is translation vector for best superposition  (result)
c  ier  - 0: a unique optimal superposition has been determined(result)
c       -1: superposition is not unique but optimal
c       -2: no result obtained because of negative weights w
c           or all weights equal to zero.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine u3b(w, x, y, n, mode, rms, u, t, ier)
      integer ip(9), ip2312(4), i, j, k, l, m1, m, ier, n, mode
      double precision w(1), x(3, 1), y(3, 1), u(3, 3), t(3), rms, sigma
      double precision r(3, 3), xc(3), yc(3), wc, a(3, 3), b(3, 3), e0, 
     &e(3), e1, e2, e3, d, spur, det, cof, h, g, cth, sth, sqrth, p, tol
     &, rr(6), rr1, rr2, rr3, rr4, rr5, rr6, ss(6), ss1, ss2, ss3, ss4, 
     &ss5, ss6, zero, one, two, three, sqrt3
      equivalence (rr(1), rr1), (rr(2), rr2), (rr(3), rr3), (rr(4), rr4)
     &, (rr(5), rr5), (rr(6), rr6), (ss(1), ss1), (ss(2), ss2), (ss(3), 
     &ss3), (ss(4), ss4), (ss(5), ss5), (ss(6), ss6), (e(1), e1), (e(2)
     &, e2), (e(3), e3)
      data sqrt3 / 1.73205080756888d+00 /
      data tol / 1.0d-2 /
      data zero / 0.0d+00 /
      data one / 1.0d+00 /
      data two / 2.0d+00 /
      data three / 3.0d+00 /
      data ip / 1, 2, 4, 2, 3, 5, 4, 5, 6 /
      data ip2312 / 2, 3, 1, 2 /
c 156 "rms.for"
      wc = zero
      rms = 0.0
      e0 = zero
      do 1 i = 1, 3
      xc(i) = zero
      yc(i) = zero
      t(i) = 0.0
      do 1 j = 1, 3
      d = zero
      if (i .eq. j) d = one
      u(i,j) = d
      a(i,j) = d
    1 r(i,j) = zero
      ier = -1
c**** DETERMINE CENTROIDS OF BOTH VECTOR SETS X AND Y
c 170 "rms.for"
      if (n .lt. 1) return 
c 172 "rms.for"
      ier = -2
      do 2 m = 1, n
      if (w(m) .lt. 0.0) return 
      wc = wc + w(m)
      do 2 i = 1, 3
      xc(i) = xc(i) + (w(m) * x(i,m))
    2 yc(i) = yc(i) + (w(m) * y(i,m))
      if (wc .le. zero) return 
      do 3 i = 1, 3
      xc(i) = xc(i) / wc
c**** DETERMINE CORRELATION MATRIX R BETWEEN VECTOR SETS Y AND X
c 182 "rms.for"
    3 yc(i) = yc(i) / wc
c 184 "rms.for"
      do 4 m = 1, n
      do 4 i = 1, 3
      e0 = e0 + (w(m) * (((x(i,m) - xc(i)) ** 2) + ((y(i,m) - yc(i)) ** 
     &2)))
c 187 "rms.for"
      d = w(m) * (y(i,m) - yc(i))
      do 4 j = 1, 3
c**** CALCULATE DETERMINANT OF R(I,J)
c 189 "rms.for"
    4 r(i,j) = r(i,j) + (d * (x(j,m) - xc(j)))
c 191 "rms.for"
      det = ((r(1,1) * ((r(2,2) * r(3,3)) - (r(2,3) * r(3,2)))) - (r(1,2
     &) * ((r(2,1) * r(3,3)) - (r(2,3) * r(3,1))))) + (r(1,3) * ((r(2,1)
     & * r(3,2)) - (r(2,2) * r(3,1))))
c**** FORM UPPER TRIANGLE OF TRANSPOSED(R)*R
c 194 "rms.for"
      sigma = det
c 196 "rms.for"
      m = 0
      do 5 j = 1, 3
      do 5 i = 1, j
      m = m + 1
c***************** EIGENVALUES *****************************************
c**** FORM CHARACTERISTIC CUBIC  X**3-3*SPUR*X**2+3*COF*X-DET=0
c 200 "rms.for"
    5 rr(m) = ((r(1,i) * r(1,j)) + (r(2,i) * r(2,j))) + (r(3,i) * r(3,j)
     &)
c 203 "rms.for"
      spur = ((rr1 + rr3) + rr6) / three
      cof = ((((((rr3 * rr6) - (rr5 * rr5)) + (rr1 * rr6)) - (rr4 * rr4)
     &) + (rr1 * rr3)) - (rr2 * rr2)) / three
c 205 "rms.for"
      det = det * det
      do 6 i = 1, 3
    6 e(i) = spur
c**** REDUCE CUBIC TO STANDARD FORM Y**3-3HY+2G=0 BY PUTTING X=Y+SPUR
c 208 "rms.for"
      if (spur .le. zero) goto 40
c 210 "rms.for"
      d = spur * spur
      h = d - cof
c**** SOLVE CUBIC. ROOTS ARE E1,E2,E3 IN DECREASING ORDER
c 212 "rms.for"
      g = (((spur * cof) - det) / two) - (spur * h)
c 214 "rms.for"
      if (h .le. zero) goto 8
      sqrth = dsqrt(h)
      d = ((h * h) * h) - (g * g)
      if (d .lt. zero) d = zero
      d = datan2(dsqrt(d),- g) / three
      cth = sqrth * dcos(d)
      sth = (sqrth * sqrt3) * dsin(d)
      e1 = (spur + cth) + cth
      e2 = (spur - cth) + sth
      e3 = (spur - cth) - sth
c.....HANDLE SPECIAL CASE OF 3 IDENTICAL ROOTS
c 224 "rms.for"
      if (mode) 10, 50, 10
c**************** EIGENVECTORS *****************************************
c 226 "rms.for"
    8 if (mode) 30, 50, 30
c 228 "rms.for"
   10 do 15 l = 1, 3, 2
      d = e(l)
      ss1 = ((d - rr3) * (d - rr6)) - (rr5 * rr5)
      ss2 = ((d - rr6) * rr2) + (rr4 * rr5)
      ss3 = ((d - rr1) * (d - rr6)) - (rr4 * rr4)
      ss4 = ((d - rr3) * rr4) + (rr2 * rr5)
      ss5 = ((d - rr1) * rr5) + (rr2 * rr4)
      ss6 = ((d - rr1) * (d - rr3)) - (rr2 * rr2)
      j = 1
      if (dabs(ss1) .ge. dabs(ss3)) goto 12
      j = 2
      if (dabs(ss3) .ge. dabs(ss6)) goto 13
   11 j = 3
      goto 13
   12 if (dabs(ss1) .lt. dabs(ss6)) goto 11
   13 d = zero
      j = 3 * (j - 1)
      do 14 i = 1, 3
      k = ip(i + j)
      a(i,l) = ss(k)
   14 d = d + (ss(k) * ss(k))
      if (d .gt. zero) d = one / dsqrt(d)
      do 15 i = 1, 3
   15 a(i,l) = a(i,l) * d
      d = ((a(1,1) * a(1,3)) + (a(2,1) * a(2,3))) + (a(3,1) * a(3,3))
      m1 = 3
      m = 1
      if ((e1 - e2) .gt. (e2 - e3)) goto 16
      m1 = 1
      m = 3
   16 p = zero
      do 17 i = 1, 3
      a(i,m1) = a(i,m1) - (d * a(i,m))
   17 p = p + (a(i,m1) ** 2)
      if (p .le. tol) goto 19
      p = one / dsqrt(p)
      do 18 i = 1, 3
   18 a(i,m1) = a(i,m1) * p
      goto 21
   19 p = one
      do 20 i = 1, 3
      if (p .lt. dabs(a(i,m))) goto 20
      p = dabs(a(i,m))
      j = i
   20 continue
      k = ip2312(j)
      l = ip2312(j + 1)
      p = dsqrt((a(k,m) ** 2) + (a(l,m) ** 2))
      if (p .le. tol) goto 40
      a(j,m1) = zero
      a(k,m1) = - (a(l,m) / p)
      a(l,m1) = a(k,m) / p
   21 a(1,2) = (a(2,3) * a(3,1)) - (a(2,1) * a(3,3))
      a(2,2) = (a(3,3) * a(1,1)) - (a(3,1) * a(1,3))
c****************** ROTATION MATRIX ************************************
c 282 "rms.for"
      a(3,2) = (a(1,3) * a(2,1)) - (a(1,1) * a(2,3))
c 284 "rms.for"
   30 do 32 l = 1, 2
      d = zero
      do 31 i = 1, 3
      b(i,l) = ((r(i,1) * a(1,l)) + (r(i,2) * a(2,l))) + (r(i,3) * a(3,l
     &))
c 288 "rms.for"
   31 d = d + (b(i,l) ** 2)
      if (d .gt. zero) d = one / dsqrt(d)
      do 32 i = 1, 3
   32 b(i,l) = b(i,l) * d
      d = ((b(1,1) * b(1,2)) + (b(2,1) * b(2,2))) + (b(3,1) * b(3,2))
      p = zero
      do 33 i = 1, 3
      b(i,2) = b(i,2) - (d * b(i,1))
   33 p = p + (b(i,2) ** 2)
      if (p .le. tol) goto 35
      p = one / dsqrt(p)
      do 34 i = 1, 3
   34 b(i,2) = b(i,2) * p
      goto 37
   35 p = one
      do 36 i = 1, 3
      if (p .lt. dabs(b(i,1))) goto 36
      p = dabs(b(i,1))
      j = i
   36 continue
      k = ip2312(j)
      l = ip2312(j + 1)
      p = dsqrt((b(k,1) ** 2) + (b(l,1) ** 2))
      if (p .le. tol) goto 40
      b(j,2) = zero
      b(k,2) = - (b(l,1) / p)
      b(l,2) = b(k,1) / p
   37 b(1,3) = (b(2,1) * b(3,2)) - (b(2,2) * b(3,1))
      b(2,3) = (b(3,1) * b(1,2)) - (b(3,2) * b(1,1))
      b(3,3) = (b(1,1) * b(2,2)) - (b(1,2) * b(2,1))
      do 39 i = 1, 3
      do 39 j = 1, 3
c****************** TRANSLATION VECTOR *********************************
c 320 "rms.for"
   39 u(i,j) = ((b(i,1) * a(j,1)) + (b(i,2) * a(j,2))) + (b(i,3) * a(j,3
     &))
   40 do 41 i = 1, 3
c********************** RMS ERROR **************************************
c 323 "rms.for"
   41 t(i) = ((yc(i) - (u(i,1) * xc(1))) - (u(i,2) * xc(2))) - (u(i,3)
     & * xc(3))
   50 do 51 i = 1, 3
      if (e(i) .lt. zero) e(i) = zero
   51 e(i) = dsqrt(e(i))
      ier = 0
      if (e2 .le. (e1 * 1.0d-05)) ier = -1
      d = e3
      if (sigma .ge. 0.0) goto 52
      d = - d
      if ((e2 - e3) .le. (e1 * 1.0d-05)) ier = -1
   52 d = (d + e2) + e1
      rms = (e0 - d) - d
      if (rms .lt. 0.0) rms = 0.0
      return 
c.....END U3B...........................................................
c----------------------------------------------------------
c                       THE END
c----------------------------------------------------------
c 338 "rms.for"
      end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Read contact information from pre-calculated files               c
c     The format of file is                                            c
c     Res  ResidueIndex  NumContacts  Contact Residues                 c
c     e.g.
c     RES    74   3   49   50   51                                     c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine readContact(contfile1,contfile2)
      PARAMETER (nmax=3000)  !Mu: maxium number of residue per protein
      PARAMETER (cmax=100)   !Mu: maximum number of contact per residue

ccc   Contact stuff added by Mu
      integer cont1(cmax,nmax),cont2(cmax,nmax)
      common/contlst/cont1,cont2  ! contact list for protein 1,2.
      common/contnum/ncont1(nmax),ncont2(nmax)
      common/contact/xcontcf  !Mu: contact cutoff

      character*100 contfile1,contfile2
      character*100 is,idu
      character*3 iaa(-1:20),iaanam,iss1(nmax),iss2(nmax),iss3(nmax)
      character*3 iss4(nmax)

****** Initializing contact arrays ********
      do i=1,nmax
         ncont1(i)=0
         ncont2(i)=0
         do j=1,cmax
            cont1(j,i)=0
            cont2(j,i)=0
         enddo
      enddo

****** read contact list files ********
      call readContFile(contfile1,ncont1,cont1)
      call readContFile(contfile2,ncont2,cont2)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Read a contact file                                              c
c     The format of file is                                            c
c     Res  ResidueIndex  NumContacts  Contact Residues                 c
c     e.g.
c     RES    74   3   49   50   51                                     c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine readContFile(contfile,ncont,cont)
      PARAMETER (nmax=3000)  !Mu: maxium number of residue per protein
      PARAMETER (cmax=100)   !Mu: maximum number of contact per residue

      character*100 contfile
      character*1000 line
      character*5,head
      integer ncont(nmax),cont(cmax,nmax)

      open(unit=15,file=contfile,status='old')
      do while (.true.)
         read(15,1001,end=100) line
         if(line(1:5).eq.'RES  ') then
            read(line,*) head, ires, ncont(ires),
     &           (cont(i,ires),i=1,ncont(ires))
c            write(*,*) head, ires, ncont(ires),
c     &           (cont(i,ires),i=1,ncont(ires))
         endif
      enddo
 100  continue
 1001 format(A1000)

      close(15)
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Read a match list file                                           c
c     The format of file is                                            c
c     Index Ch1 Resid1 AA1 Ch2 Resid2                                  c 
c     e.g.                                                             c
c     1 C    31 THR C    31 THR                                        c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine readMatchFile(matchfile,matresch1, matresid1,
     &     matresch2,matresid2,nmatch)
      implicit none
      integer nmax
      PARAMETER (nmax=3000)  !Mu: maxium number of residue per protein


      character*500 matchfile
      character*100 line
      character*3   AA1,AA2
      character*5   matresid1(nmax),matresid2(nmax)
      character     matresch1(nmax),matresch2(nmax)

      integer i,ires,nmatch

      open(unit=16,file=matchfile,status='old')
      i=0
 5    read(16,1001,end=10)line

      if(line(1:1).eq.'#')     goto 5
      if(line(1:5).eq.'Index') goto 5

      i=i+1
      read(line,'(I5)')ires

      if(ires.eq.i)then
         read(line,205)ires,matresch1(i),matresid1(i),AA1,
     &        matresch2(i),matresid2(i),AA2
c         write(*,*)ires, matresch1(i), matresid1(i),AA1,
c     &        matresch2(i),matresid2(i),AA2
         nmatch=i
      endif
      goto 5
 10   continue

 205  format(i5,1x,a1,1x,a5,1x,a3,1x,a1,1x,a5,1x,a3)
 1001 format(A100)

      close(16)
      return
      end

****************************************************************
*     initialize alignment mapping  arrays
*     the array map original position to aligned position
****************************************************************
      subroutine init_ali(ali1,ali2,nseq1,nseq2)
      parameter(nmax=3000)

      integer ali1(nmax), ali2(nmax)
      integer nseq1,nseq2,i

ccc   initialize alignment mapping array
      do i=1,nseq1
         ali1(i)=-1
      enddo
      do i=1,nseq2
         ali2(i)=-2
      enddo

      return
      end





****************************************************************
*     calculate contact overlap factor
****************************************************************
      subroutine get_fcol(pos1,pos2,ali1,ali2,imap1,imap2,fcol,col)
      parameter(nmax=3000)
      parameter(cmax=100)

      integer imap1(nmax),imap2(nmax)  !Mu: Map aligned position to original position
      integer ali1(nmax), ali2(nmax)   !Mu: Map original poistion to aligned position

ccc   Contact stuff added by Mu
      integer cont1(cmax,nmax),cont2(cmax,nmax)
      common/contlst/cont1,cont2
      common/contnum/ncont1(nmax),ncont2(nmax)
      integer pos1,pos2   !Mu: position in alignment, typically same numbers
      integer col
      real fcol

      col=0
      do 31, kk1=1,ncont1(imap1(pos1))
         ii=cont1(kk1,imap1(pos1))
         !print *, 'contact list prot1: ',imap1(pos1),'-',ii
         if(ali1(ii).gt.0)then ! ii must be in aligned region
            do 41, kk2=1,ncont2(imap2(pos2))
               jj=cont2(kk2,imap2(pos2))
               !print *, 'contact list prot2: ',imap2(pos2),'-',jj
               if(ali2(jj).eq.ali1(ii))then
                  col=col+1
                  !print *, imap1(pos1),'-',ii,imap2(pos2),'-',jj
                  goto 31
               endif
 41         continue
         endif
 31   continue

      fcol=0.0
      if(ncont1(imap1(pos1)).gt.0.and.ncont2(imap2(pos2)).gt.0)then
         fcol=(real(col)/ncont1(imap1(pos1))+real(col)/ncont2(imap2(pos2)))/2
      endif

      return
      end

cccccc================================================================
      subroutine get_fcol_4s_old(pos1,pos2,is_ali,fcol,col)
      parameter(nmax=3000)
      parameter(cmax=100)

      integer is_ali(nmax, cmax)   !alignment map

ccc   Contact stuff added by Mu
      integer cont1(cmax,nmax),cont2(cmax,nmax)
      common/contlst/cont1,cont2
      common/contnum/ncont1(nmax),ncont2(nmax)
      integer pos1,pos2   !position on sequence 1 and 2
      real fcol
      integer col

      col=0
      fcol=0
      if(is_ali(pos1,pos2).lt.0) return !do not continue if i,j are not aligned
      do 31, kk1=1,ncont1(pos1)
         ii=cont1(kk1,pos1)
         do 41, kk2=1,ncont2(pos2)
            jj=cont2(kk2,pos2)
            if(is_ali(ii,jj).gt.0)then
               col=col+1
            endif
 41      continue
 31   continue


      if(ncont1(pos1).gt.0.and.ncont2(pos2).gt.0)then
         fcol=real(col)/real(ncont1(pos1))
         fcol=(fcol+real(col)/real(ncont2(pos2)))/2
      endif

      return
      end

ccccc================================================ccccc
ccccc generate the alignment in sequences for output ccccc
ccccc================================================ccccc
      subroutine get_alnseq(alnmap1,alnmap2,laln,
     &     seq1,seq2,nchain1,nchain2,cterm1,cterm2,
     &     alnx1,alny1,alnz1,alnx2,alny2,alnz2,d_cutoff,
     &     alnseq1,alnseq2,markseq,lseq)
      parameter(maxr=3000)
      parameter(maxk=20)

      character seq1(maxr),seq2(maxr),aa
      character fseq1(maxr),fseq2(maxr)
      character alnseq1(maxr),alnseq2(maxr)
      character markseq(maxr)

      real alnx1(maxr),alny1(maxr),alnz1(maxr)
      real alnx2(maxr),alny2(maxr),alnz2(maxr)
      real d_cutoff  ! distance cutoff for marks

      integer alnmap1(maxr),alnmap2(maxr)
      integer cterm1(0:maxk),cterm2(0:maxk)
      integer laln,lseq
      integer pos    !indexes on final aligned sequence


      lseq=cterm1(nchain1)+cterm2(nchain2)-laln

      do i=1,lseq
         alnseq1(i) = '-'
         alnseq2(i) = '-'
         markseq(i) = ' '
      enddo

      call formatSeq(nchain1,cterm1,seq1,fseq1)
      call formatSeq(nchain2,cterm2,seq2,fseq2)

      pos=0     
      i1_prev=0  !position in unaligned sequence
      i2_prev=0

      alnmap1(laln+1)=cterm1(nchain1)+1  ! artificial alignment for handling sequence ends properly
      alnmap2(laln+1)=cterm2(nchain2)+1
      do 200 k=1,laln+1
         i1=alnmap1(k)
         i2=alnmap2(k)

         do j=i1_prev+1,i1-1
            pos=pos+1
            alnseq1(pos)=fseq1(j)
         enddo

         do j=i2_prev+1,i2-1
            pos=pos+1
            alnseq2(pos)=fseq2(j)
         enddo

         if(k.eq.laln+1)goto 200  !ignore the last artificial alignment

         pos=pos+1
         alnseq1(pos)=fseq1(i1)
         alnseq2(pos)=fseq2(i2)

         dist=sqrt( (alnx1(k)-alnx2(k))**2
     &        +     (alny1(k)-alny2(k))**2
     &        +     (alnz1(k)-alnz2(k))**2 )

         if(dist.le.d_cutoff)then
            markseq(pos)=':'
         endif

         i1_prev=i1
         i2_prev=i2
 200  continue

      return
      end


ccccc==================================================ccccc
ccccc generate sequential order independent alignment  ccccc
ccccc==================================================ccccc
      subroutine get_matchlst(nseq1,nseq2,mf1,mf2,xtmf1,xtmf2,
     &     ytmf1,ytmf2,ztmf1,ztmf2,respdbid1,respdbid2,seq1,seq2,
     &     ncont1,ncont2,chname1,chname2,naln,d_output,matchlst)

      implicit none
      integer maxr,maxk
      parameter(maxr=3000)
      parameter(maxk=20)

      common/d0/d0,anseq

      character*5 respdbid1(maxr),respdbid2(maxr)
      character*3 seq1(maxr),seq2(maxr)
      character   chname1(maxr), chname2(maxr)
      character*80 matchlst(maxr)
      character*3  note

      integer mf1(maxr),mf2(maxr)
      real xtmf1(maxr),ytmf1(maxr),ztmf1(maxr)
      real xtmf2(maxr),ytmf2(maxr),ztmf2(maxr)
      real d0,anseq,d_output

      integer id_chain1(maxr),id_chain2(maxr)  !original index to chain ID
      integer ncont1(maxr),ncont2(maxr)      !number of contacts
      integer nseq1, nseq2, naln
      integer is_ali(maxr, maxr)   !alignment map
      integer col,ncol
      real    fcol,tms,tms1,ss,dis
      integer i,j,k,i1,i2,k1,k2

      call init_ali_map( is_ali, nseq1, nseq2 )
      do k=1,naln
         i=mf1(k)
         j=mf2(k)
         is_ali(i,j) = 1
      enddo

      tms=0
      tms1=0
      do 100 k=1,nseq2

         i=-1
         do j=1,naln
            i1=mf1(j)
            i2=mf2(j)
            if(i2.eq.k)then
               i=j
               goto 10
            endif
         enddo
 10      continue

         if(i.eq.-1)then
            write(matchlst(k),204)k,chname2(k),
     &           respdbid2(k),seq2(k),' -'
            goto 100
         endif

         note='   '
         dis=sqrt((xtmf1(i)-xtmf2(i))**2+
     &        (ytmf1(i)-ytmf2(i))**2+(ztmf1(i)-ztmf2(i))**2)
         if(dis.lt.d_output) note(2:2) = ':'

         call get_fcol_4s_old(i1,i2,is_ali,fcol,col)
         if(seq1(i1).eq.seq2(i2))then
            note(3:3) = '*'
         endif

         write(matchlst(k),205)k,chname1(i1),
     &        respdbid1(i1),seq1(i1),chname2(i2),
     &        respdbid2(i2),seq2(i2),dis,col,
     &        ncont1(i1),ncont2(i2),note

         ss=1/(1+(dis/d0)**2)
         tms=tms+ss
         tms1=tms1+fcol*ss
         !print *,'fcol=',fcol,' ss=',ss

 100  continue
      tms=tms/anseq
      tms1=tms1/anseq
!      write(*,'(a,f5.3,a,f5.3)')'TMS = ',tms,
!     &     ', TMS (cont) = ',tms1

 204  format(i5,20x,a3,1x,a5,2x,a3,4x,21x,2a)
 205  format(i5,2x,2(a3,1x,a5,2x,a3,4x),f7.3,i5,i4,i4,a4)

      return
      end


ccccc=========================================================ccccc
ccccc Convert one-letter amino acid code to three-letter code ccccc
ccccc=========================================================ccccc
      subroutine AA1toAA3(resnm3,resnm1)
      implicit none

      character*3 aa(-1:19),resnm3
      data aa/ 'UNK','GLY','ALA','SER','CYS',
     &     'VAL','THR','ILE','PRO','MET',
     &     'ASP','ASN','LEU','LYS','GLU',
     &     'GLN','ARG','HIS','PHE','TYR',
     &     'TRP'/

      character*1 slc(-1:19),resnm1
      data slc/'X','G','A','S','C',
     &     'V','T','I','P','M',
     &     'D','N','L','K','E',
     &     'Q','R','H','F','Y',
     &     'W'/

      integer j

      resnm3='UNK'
      do j=-1,19
         if(resnm1.eq.slc(j))then
            resnm3=aa(j)
            return
         endif
      enddo

      return
      end


ccccc=========================================================ccccc
ccccc Format sequence: odd chains are lower case, even upper  ccccc
ccccc=========================================================ccccc

      subroutine formatSeq(nchain,cterm,seq,fseq)
      integer maxr,maxk
      parameter(maxr=3000)
      parameter(maxk=20)

      character*1 seq(maxr),fseq(maxr),aa
      integer cterm(0:maxk)
      integer nchain,ich,i,j

      ! make a copy of the original sequence
      do i=1,maxr
         fseq(i)=seq(i)
      enddo

      ! first make the sequence all upper case
      do i=1,nchain
         do j=cterm(i-1)+1, cterm(i)
            ich=ichar(seq(j))
            if(ich.gt.ichar('Z'))then
               aa =char(ich-32)
               fseq(j)=aa
            endif
         enddo
      enddo

      ! then make sequence of even chains lower case
      do i=1,nchain
         do j=cterm(i-1)+1, cterm(i)
            if(mod(i,2).eq.0)then
               ich=ichar(seq(j))
               aa =char(ich+32)
               fseq(j)=aa
            endif
         enddo
      enddo

      return
      end
ccccc==========================================================cccc


      subroutine getBasename(fullname,basename)
      implicit none
      character fullname*(*),basename*(*)
      integer i,j,l,n_sta,n_end

      l = len(fullname)
      n_sta=0
      n_end=l
      do i=1,l
         if(fullname(i:i).eq.'/')then
            n_sta=i
         else if(fullname(i:i).ne.' ')then
            n_end=i
         endif
      enddo

      l = n_end - n_sta
      do i=1,l
         j=n_sta+i
         basename(i:i)=fullname(j:j)
      enddo

      ! maximum length of basename is 20 characters
      do i=l+1,25
         basename(i:i)=' '
      enddo

      return
      end

cccccc==================================================
ccccc initialize a matrix
ccccc===================================================


      subroutine init_ali_map( is_ali, nseq1, nseq2 )
      parameter(maxr=3000)

      integer nseq1,nseq2
      integer is_ali(maxr, maxr)   !alignment map
      integer i,j

ccc   is_ali(i,j) indicate whether i, j are in alignment dis<d_col
      do i=1,nseq1
         do j=1,nseq2
            is_ali(i,j) = -1
         enddo
      enddo

      return
      end

ccccc===================================================ccccc
ccccc  P-value calculation for IS-score.                ccccc
ccccc  Parameters were obtained empirically from        ccccc
ccccc  random interface alignments                      ccccc
ccccc===================================================ccccc
      subroutine calcPvalue(score,len,pvalue,z)
      implicit none

      real*8  loc, scale, logl
      real*8  z, pvalue
      real    score
      integer len,len1

      len1=len
      if(len1.lt.20)len1=20
      logl = log(dble(len1))

      pvalue=-1
      z=-1

      if( len1 .lt. 55 ) then
         loc   = 0.0806 + 0.0034*logl
         scale = 0.0277 - 0.0040*logl
      else
         loc   = 0.0794 + 0.0033*logl
         scale = 0.0329 - 0.0054*logl
      endif

      if(scale.lt.0.002) scale=0.002 !prevent insane values

      call calcEVDPV(score,loc,scale,pvalue,z)


      return     
      end 



ccccc===================================================ccccc
ccccc  P-value calculation for Gumbel Distribution.     ccccc
ccccc===================================================ccccc
      subroutine calcEVDPV(score,loc,scale,pvalue,z)

      real*8  loc, scale
      real*8  z, pvalue
      real    score

      z = (score - loc) / scale
      if( z .lt. 35 ) then
         pvalue = 1 -  dexp( - dexp( -z ) ) !Gumbel distribution
      else
c        below the double precision limit, use taylor expansion approaximation
         pvalue = dexp( -z )
      endif

      return     
      end 


ccccc===================================================ccccc
ccccc  Print the help message                           ccccc
ccccc===================================================ccccc
      subroutine printHelp()

         write(*,*)
         write(*,'(A)')'Usage: IS-score <options> pdbfile1 pdbfile2 '//
     &        'contactfile1 contactfile2'
         write(*,*)
         write(*,'(A)')'By default, pdbfile1 is a model (or template), '//
     &        'and pdbfile2 is the target of the model.'
         write(*,'(A)')'The length of the latter is'//
     &         ' used for score normalization.'

         write(*,*)
         write(*,'(A)')'Options:'
         write(*,'(A)')'     -t        Use TM-score as the similarity'//
     &        ' measure. The default is IS-score.'

         write(*,'(A)')'     -eqpos    Residues at the same position '//
     &        'is considered equivilent regardless of their identities.'
         write(*,'(A)')'     -d <num>  Normalize the score with a fixed '//
     &        'distance factor d0.'

         write(*,'(A)')'     -v        0 - no alignment, '//
     &     ' 1 - concise alignment, 2 - detailed alignment.'
         write(*,*)

      return
      end
