       integer ncol,ntm,i,j,k,l,t,x,rowno,count,LIST,maxrow
       integer minrow,col,avgrow,rowflag,redunt,num_times
       DOUBLE PRECISION::tol,tm
       DOUBLE PRECISION, ALLOCATABLE,DIMENSION(:,:,:)::CORR,bionet
       DOUBLE PRECISION, ALLOCATABLE,DIMENSION(:,:,:)::ALIGN
       DOUBLE PRECISION, ALLOCATABLE,DIMENSION(:,:) ::C1C2avg,AVG
       DOUBLE PRECISION, ALLOCATABLE,DIMENSION(:) ::DIFF
       DOUBLE PRECISION, ALLOCATABLE,DIMENSION(:,:,:,:) ::covar
       integer, ALLOCATABLE,DIMENSION(:) ::AVGDIFF,nrow
       character (len = 12),ALLOCATABLE, dimension(:) :: file_i
       character (len = 12) :: filename
cc       CHARACTER (len=100) :: t_string
       character(len=5),ALLOCATABLE,DIMENSION(:) :: t_i
cc       integer :: pos1 = 1, pos2, n

cc       call system("wc -m <time.txt > t_len")
       call system("wc -l <time.txt > num_times")
cc       open(unit=11,file='t_len',status='old')
cc       read(11,*)t_len
cc       close(11)
       open(unit=11,file='num_times',status='old')
        read(11,*)num_times
       close(11)


       ALLOCATE(t_i(num_times))
       ALLOCATE(file_i(num_times+1))
cc       ALLOCATE(t_string(t_len))
       open(unit=11,file='time.txt',status='old',action='read')
       do j=1,num_times
          read(11,*)t_i(j)
       end do
cc       read(unit=11,fmt='(A)')t_string
       close(11)

       
cc This segment of code below was taken from rosettacode.org.  The section is "Tokenize a string", subsection "Fortran"
cc       DO n = 1, num_times
cc          pos2 = INDEX(t_string(pos1:), " ")
cc          IF (pos2 == 0) THEN
cc             t_i(n) = t_string(pos1:)
cc             EXIT
cc          END IF
cc          t_i(n) = t_string(pos1:pos1+pos2-2)
cc          pos1 = pos2+pos1
cc       END DO

c Opening the file "avg_NKG2D_NKLs_new.txt" and assigning it's contents to "inpl"
       call system("wc -l <avg_new.txt > inpl")
c Asign something to "inpc" using something from avg_NKG2D...txt
       call system("awk '(FNR==2){print NF}' 
     &      avg_new.txt> inpc")
   

c Read file "inpl", make it "avgrow" variable
       open(unit=11,file='inpl',status='old')
       read(11,*)avgrow
       close(11)
       
c Read file "inpc", make it "ncol" variable
       open(unit=11,file='inpc',status='old')
       read(11,*)ncol
       close(11)

c Allocate?  Not sure what this is.  Some sort of assigning a new matrix?
       ALLOCATE(AVG(avgrow,ncol))
       ALLOCATE(AVGDIFF(avgrow-1))
       ALLOCATE(nrow(avgrow))
       AVG=0.0d0
       AVGDIFF=0
       
c Open avg_NKG2D_NKLs_new.txt and read it (file ID = 12)
       open(unit=12,file='avg_new.txt',
     &        status='old')
       do i=1,avgrow
        read(12,*)(AVG(i,j),j=1,ncol)
       enddo
       close(12)
  
c For every row, do this thing.  Some sort of percent difference in averages?

c Write this percent difference(?) matrix to a file
c        write(*,*)(AVGDIFF(i),i=1,avgrow-1)

c tol = tolerance?  Or something else?   
        tol=1e+8
c Remove "inpl"
        call system("rm inpl")

c Create all these files and do something where they are writing and being written to at the same time?
        do i=1, num_times+1
           if (i == (num_times+1)) then
              file_i(i) = "zzz"
              exit
           else
              write(filename,*)trim(t_i(i)), "_min.txt"
              filename=adjustl(filename)
              file_i(i) = trim(filename)
              call system("wc -l <" // filename // ">> 
     & inpl")
           end if
        enddo

       nrow=0

c Do the same thing we did earlier with 'inpl', but assign it to this varialbe 'nrow'?
       open(unit=11,file='inpl',status='old')
       do i=1,avgrow
       read(11,*)nrow(i)
       enddo
       close(11) 

c Assign maxrow/minrow variables from the max/min values from each row
        maxrow=MAXVAL(nrow)
        minrow=MINVAL(nrow)
        do i=1,avgrow
          if(minrow.eq.nrow(i)) then
            rowflag=i
          endif
        enddo
c Print "smallest set = [Row with minrow]"
        write(*,*)"smallest set =",rowflag

c Allocate again...?  Also don't know if these are built-in functions (bionet, ALIGN, etc.)
         ALLOCATE(bionet(avgrow,maxrow,ncol))
         ALLOCATE(ALIGN(avgrow,minrow,ncol))
         ALLOCATE(C1C2avg(avgrow,ncol))
         ALLOCATE(CORR(avgrow,ncol,ncol))
         ALLOCATE(covar(ncol,ncol,avgrow,avgrow))

         bionet=0.0d0
         ALIGN=0.0d0
         C1C2avg=0.0d0
         CORR=0.0d0
         covar=0.0d0
         LIST=0

c Open and read these files (to bionet?) using file ID 13 
         x=1
         do 
            if (file_i(x) == "zzz") exit
         open(unit=13,file=file_i(x),
     &                                   status='old')
         do i=1,nrow(x)
         read(13,*)(bionet(x,i,j),j=1,ncol)
         enddo
         close(13)
         x=x+1
         enddo

c       open(unit=13,file='IL2_NKG2D_128min_live_NKLs_spaces_sayak.txt',
c     &                                   status='old')
c         do i=1,nrow(8)
c         read(13,*)(bionet(8,i,j),j=1,ncol)
c         enddo
c         close(13)

c Revert values that are less than 0 to 0
        do k=1,avgrow
          do j=1,ncol
           do i=1,nrow(k)
           if(bionet(k,i,j).lt.0) then
             bionet(k,i,j)=0.0d0
            endif
           enddo
          enddo
         enddo

c Don't know what this is doing.  Calculating covariance?
         do k=1,avgrow
          do i=1,ncol
           do rowno=1,nrow(k)
            C1C2avg(k,i)=C1C2avg(k,i)+bionet(k,rowno,i)
           enddo
          enddo
         enddo

         do k=1,avgrow
          do j=1,ncol
            C1C2avg(k,j)=C1C2avg(k,j)/(1.0*nrow(k))
          enddo
         enddo
   
c Continuing covariance calculations?  Maybe 'bionet' is a variance?
         do k=1,avgrow
          do i=1,ncol
           do j=1,ncol
            do rowno=1,nrow(k)
            CORR(k,i,j)=CORR(k,i,j)+bionet(k,rowno,i)*bionet(k,rowno,j)
            enddo
           enddo
          enddo
         enddo

         do k=1,avgrow
          do i=1,ncol
           do j=1,ncol
            CORR(k,i,j)=CORR(k,i,j)/(1.0*nrow(k))
           enddo
          enddo
         enddo

c Continuing calculations (I know they shouldn't be terribly important to me to get the simulation to work, i.e. I shouldn't have to change anything)
         do k=1,avgrow
           do i=1,ncol
            do j=1,ncol
           CORR(k,i,j)=CORR(k,i,j)-C1C2avg(k,i)*C1C2avg(k,j)
           enddo
          enddo
         enddo

c Open a new file and write the correlation coefficients to it (file ID=60)
         open(unit=60,status='unknown')
         do k=1,avgrow
           do i=1,ncol
             do j=1,ncol
              
              write(60,'(10F32.16)')CORR(k,i,j)
             
            enddo
           enddo
         enddo
         close(60)

c Open a new file and write the covariance(?) to it (file ID=61)
          open(unit=61,status='unknown')
              do k=1,avgrow
                write(61,'(512F32.16)')(C1C2avg(k,i),i=1,ncol)
              enddo
           close(61)

           end

