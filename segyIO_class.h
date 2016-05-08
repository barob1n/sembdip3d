int  checkEndian(void);
void headerAscii(unsigned char *s1, unsigned char *s2, int n);
int  segyWriteHeader(FILE *f, char *txt_hdr, bhed *bh, int endian);
int  segyWriteTrace(FILE *f, segy *tr_hdr, float *tr_data,
	int ns, int endian);
void swap4bytes(char *b);
void swap2bytes(char *b);
int  getHeaderSwap(char *p1, int byte, int count, int endian);
void swapBhead(bhed *bh);
//void swapTrhead(char *bchr);
void swapTrHead (segy *bchar);
int  segyReadHeader(FILE *f, char *txt_hdr, bhed *bh, int endian);
int  segyReadTrace(FILE *f, bhed *bh, segy *tr_hdr, float *tr_data,
	int ns, int endian);
void toebc(unsigned char *s1, unsigned char *s2, int n);
void makeTapeHdr(unsigned char *s1);
void makeBinHdr(bhed *bhdr, int ns, float dt);
void makeTrHdr(segy *thdr, int ns, float dt, int scal);
void printBinHdr(bhed *bhdr);
