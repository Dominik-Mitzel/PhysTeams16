
 void createNtuple(){
   gROOT->Reset();

   FILE *fp = fopen("debug.dat","r");

   Float_t w,z1,z2,z3,z4,sqrts, sqrtt, pt;
   Int_t ncols;
   Int_t nlines = 0;
   TFile *f = new TFile("basicTuple.root","RECREATE");
   TNtuple *ntuple = new TNtuple("ntuple","data from file","w:z1:z2:z3:z4:sqrts:sqrtt:pt");

   while (1) {
     ncols = fscanf(fp,"%f %f %f %f %f %f %f %f",&w,&z1,&z2,&z3,&z4,&sqrts,&sqrtt,&pt);
     if (ncols < 0) break;    
     //if (nlines < 5) printf("var1=%8f, var2=%8f\n",var1,var2);
     //h1->Fill(var1);
     ntuple->Fill(w,z1,z2,z3,z4,sqrts,sqrtt,pt);
     nlines++;
   }
   printf(" found %d points\n",nlines);
   
   fclose(fp);

   f->Write();
}

  
  
