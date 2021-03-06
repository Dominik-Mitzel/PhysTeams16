
 void createNtuple(){
   gROOT->Reset();

   FILE *fp = fopen("kinematics.dat","r");

   Float_t w,z1,z2,z3,z4,sqrts, sqrtt, pt, x1, x2;
   Int_t ncols;
   Int_t nlines = 0;
   TFile *f = new TFile("basicTuple.root","RECREATE");
   TNtuple *ntuple = new TNtuple("ntuple","data from file","w:z1:z2:z3:z4:sqrts:sqrtt:pt:x1:x2");

   while (1) {
     ncols = fscanf(fp,"%f %f %f %f %f %f %f %f %f %f",&w,&z1,&z2,&z3,&z4,&sqrts,&sqrtt,&pt,&x1,&x2);
     if (ncols < 0) break;    
     ntuple->Fill(w,z1,z2,z3,z4,sqrts,sqrtt,pt,x1,x2);
     nlines++;
   }
   printf(" found %d points\n",nlines);
   
   fclose(fp);

   f->Write();
}

  
  
