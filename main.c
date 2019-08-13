#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define LUNGIME_NUME_FISIER 32
#define HEADER_SIZE 54

typedef struct ImgData
{
    unsigned char *header;
    unsigned int W;
    unsigned int H;
    unsigned int size;
    unsigned char *data;
} ImgData;

typedef struct dreptunghi
{
    int stanga;
    int dreapta;
    int sus;
    int jos;
} dreptunghi;

typedef struct matchData
{
    int cifra;
    int x;
    int y;
    double corr;
} matchData;

void ImgDatatoBMP(char const filename[], ImgData image)
{
    FILE *f_exBMP = fopen(filename, "wb");
    if(f_exBMP == NULL)
    {
        printf("Nu se poate deschide fisierul %s\n", filename);
        exit(EXIT_FAILURE);
    }

    int padding;
    if(image.W % 4 != 0)
        padding = 4 - (3 * image.W) % 4;
    else
        padding = 0;
    unsigned char paddingarray[3] = {0};

    fwrite(image.header, sizeof(unsigned char), 54, f_exBMP); //Scrie header
    for(unsigned int i = 0; i < image.H; i++) //prima linie din vector este de fapt ultima linie de pixeli a pozei
    {
        fwrite(image.data + (image.H - i - 1)*image.W *3, sizeof(unsigned char), image.W * 3, f_exBMP);
        fwrite(paddingarray, sizeof(unsigned char), padding, f_exBMP);
    }
    fflush(f_exBMP);
    fclose(f_exBMP);
}

void BMPtoImgData(char const filename[], ImgData *image)
{
    FILE *f_BMP = fopen(filename, "rb");
    if(f_BMP == NULL)
    {
        printf("Eroare la deschiderea fisierului %s\n", filename);
        exit(EXIT_FAILURE);
    }

    image->header = (unsigned char *) malloc(sizeof(unsigned char) * HEADER_SIZE);
    if(image->header == NULL)
    {
        printf("Nu se poate aloca memorie pentru header-ul imaginii %s\n", filename);
        exit(EXIT_FAILURE);
    }

    fread(image->header, sizeof(char), 54, f_BMP);

    //Resetare + Skip primilor 18 octeti
    fseek(f_BMP, 18, SEEK_SET);

    //Width
    fread(&(image->W), sizeof(unsigned int), 1, f_BMP);

    //Height
    fread(&(image->H), sizeof(unsigned int), 1, f_BMP);

    image->size = image->H * image->W;

    //Resetare + mutare pe pozitia 54 (dupa header)
    fseek(f_BMP, 54, SEEK_SET);

    image->data = (unsigned char*) malloc(image->size * 3); // marime * 3 canale de culoare
    if(image->data == NULL)
    {
        printf("Eroare la alocarea memoriei necesare liniarizarii imaginii.\n");
        exit(EXIT_FAILURE);
    }

    int padding;
    if(image->W % 4 != 0)
        padding = 4 - (3 * image->W) % 4;
    else
        padding = 0;

    for(int i = 0; i < image->H; i++) //Citire de jos in sus si de la stanga la dreapta
    {
        //Pointer in fisier la inceputul ultimei linii
        fseek(f_BMP, (-1)*((i+1)*(image->W * 3 + padding)), SEEK_END);
        //prima linie din vector corespunde ultimei linii de pixeli din imagine
        fread(image->data + i*image->W *3, sizeof(unsigned char), image->W*3, f_BMP);
    }
    fclose(f_BMP);
}

unsigned int * Genereaza_Permutare(unsigned int *rand, unsigned int size)
{
    unsigned int r, aux, i;
    unsigned int *perm;
    perm = (unsigned int*) malloc(sizeof(unsigned int) * size);
    if(perm == NULL)
    {
        printf("Eroare la alocarea memoriei pentru array-ul de permutari!\n");
        exit(EXIT_FAILURE);
    }
    for(i = 0; i < size; i++)
        perm[i] = i;
    for(i = size-1; i >= 1; i--)
    {
        r = rand[size-i] % (i+1);
        aux = perm[r];
        perm[r] = perm[i];
        perm[i] = aux;
    }
    return perm;
}


void XORSHIFT32(unsigned int **v, unsigned int size, unsigned int seed)
{
    unsigned int x;
    unsigned int i;
    (*v)[0] = seed;
    x = seed;
    for (i = 1; i < size; i++)
    {
        x = x ^ x << 13;
        x = x ^ x >> 17;
        x = x ^ x << 5;
        (*v)[i] = x;
    }
}

void Inverseaza_Permutare(unsigned int **v, unsigned int size)
{
    int i;
    unsigned int *perm;
    perm = (unsigned int*) malloc(sizeof(unsigned int) * size);
    if(perm == NULL)
    {
        printf("Eroare la alocarea memoriei pentru Generarea permutarilor!");
        exit(EXIT_FAILURE);
    }
    for(i = 0; i < size; i++)
    {
        perm[i] = (*v)[i];
    }
    for(i = 0; i < size; i++)
    {
        (*v)[perm[i]] = i;

    }
    free(perm);
}

void GetBytes(unsigned int x, unsigned char **byte)
{
    (*byte)[2] = (x >> 16) & 255;
    (*byte)[1] = (x >> 8) & 255;
    (*byte)[0] = x & 255;
}


void PixelXOR(unsigned char **v, unsigned char *a, unsigned char *b, unsigned char *c)
{
    (*v)[0] = a[0] ^ b[0] ^ c[0];
    (*v)[1] = a[1] ^ b[1] ^ c[1];
    (*v)[2] = a[2] ^ b[2] ^ c[2];
}

void Decripteaza(char const in_filename[], char const out_filename[], char const secretCode_filname[])
{
    ImgData image;
    BMPtoImgData(in_filename, &image);

    unsigned int r0, sv;
    FILE *f_scode = fopen(secretCode_filname, "r");
    if(f_scode == NULL)
    {
        printf("Eroare la deschiderea fisierului %s\n", secretCode_filname);
        exit(EXIT_FAILURE);
    }
    fscanf(f_scode, "%u %u", &r0, &sv); //codurile secrete

    unsigned int *rand = (unsigned int*) malloc(image.size * 2 * sizeof(unsigned int));
    if(rand == NULL)
    {
        printf("Nu se poate aloca memorie pentru array-ul de numere random!");
        exit(EXIT_FAILURE);
    }
    XORSHIFT32(&rand, image.size * 2, r0); //Genereaza vector rand de numere random

    unsigned int *perm;
    perm = Genereaza_Permutare(rand, image.size);
    if(perm == NULL)
    {
        printf("Nu se poate aloca memorie pentru array-ul de permutari!\n");
    }
    Inverseaza_Permutare(&perm, image.size);

    //Decriptare

    unsigned char *SVbyte = malloc(sizeof(unsigned char) * 3);
    if(SVbyte == NULL)
    {
        printf("Nu se poate aloca memorie pentru array-ul temporar de descompunere in biti!\n");
        exit(EXIT_FAILURE);
    }

    unsigned char *RNDbyte = malloc(sizeof(unsigned char) * 3);
    if(RNDbyte == NULL)
    {
        printf("Nu se poate aloca memorie pentru array-ul temporar de descompunere in biti!\n");
        exit(EXIT_FAILURE);
    }

    unsigned char *CopieIMG = (unsigned char*) malloc(sizeof(unsigned char) * image.size * 3); //Imaginea dupa permutare
    if(CopieIMG == NULL)
    {
        printf("Nu se poate aloca memorie pentru copia imaginii!");
        exit(EXIT_FAILURE);
    }
    for(int i = 0; i < image.size * 3; i=i+3)
    {
        CopieIMG[i] = image.data[i];
        CopieIMG[i+1] = image.data[i+1];
        CopieIMG[i+2] = image.data[i+2];
    }

    //Conversie numar in vector de bytes (fara byte-ul cel mai semnificativ)
    GetBytes(sv, &SVbyte);
    GetBytes(rand[image.size], &RNDbyte);

    //Primul parametru primeste ceilalte 3 parametrii XOR-ati.
    PixelXOR(&(image.data), SVbyte, image.data, RNDbyte);
    unsigned char *aux;
    unsigned char *aux2;
    for(int index = 3; index < image.size*3; index=index+3)
    {
        GetBytes(rand[image.size+(index/3)], &RNDbyte);
        aux = CopieIMG + index;
        aux2 = image.data+index;
        PixelXOR(&aux2, aux-3, aux, RNDbyte);

    }
    //Permuta pixelii decriptati

    unsigned char *auxx = (unsigned char*) malloc(sizeof(unsigned char) * image.size * 3);
    if(auxx == NULL)
    {
        printf("Nu se poate aloca memorie pentru copia imaginii decriptate!");
        exit(EXIT_FAILURE);
    }
    for(int i = 0; i < image.size * 3; i++)
    {
        auxx[i] = image.data[i];
    }

    for(int i = 0; i < image.size * 3; i=i+3)
    {
        image.data[perm[i/3]*3] = auxx[i]; //Pixelul i va fi pe perm[i]
        image.data[perm[i/3]*3+1] = auxx[i+1];
        image.data[perm[i/3]*3+2] = auxx[i+2];
    }

    ImgDatatoBMP(out_filename, image);
    free(auxx);
    free(SVbyte);
    free(RNDbyte);
    free(CopieIMG);
    free(perm);
    free(rand);
    free(image.data);
    free(image.header);
}

void Cripteaza(char const in_filename[], char const out_filename[], char const secretCode_filname[])
{
    ImgData image;
    BMPtoImgData(in_filename, &image);
    unsigned int r0, sv;

    FILE *f_scode = fopen(secretCode_filname, "r");
    if(f_scode == NULL)
    {
        printf("Eroare la deschiderea fisierului %s\n", secretCode_filname);
        exit(EXIT_FAILURE);
    }
    fscanf(f_scode, "%u %u", &r0, &sv); //codurile secrete

    unsigned int *rand = (unsigned int*) malloc(image.size * 2 * sizeof(unsigned int));
    if(rand == NULL)
    {
        printf("Nu se poate aloca memorie pentru array-ul de numere random!");
        exit(EXIT_FAILURE);
    }
    XORSHIFT32(&rand, image.size * 2, r0); //Genereaza vector rand de numere random

    unsigned int *perm;
    perm = Genereaza_Permutare(rand, image.size);
    if(perm == NULL)
    {
        printf("Nu se poate aloca memorie pentru array-ul de permutari!\n");
        exit(EXIT_FAILURE);
    }

    //Pixelii vectorului initial se permuta conform permutarii perm[]
    unsigned char *IMGPerm = (unsigned char*) malloc(sizeof(unsigned char) * image.size * 3); //Imaginea dupa permutare
    if(IMGPerm == NULL)
    {
        printf("Nu se poate aloca memorie pentru array-ul de pixeli permutati!");
        exit(EXIT_FAILURE);
    }
    for(int i = 0; i < image.size * 3; i=i+3)
    {
        IMGPerm[perm[i/3]*3] = image.data[i]; //Pixelul i va fi pe perm[i]
        IMGPerm[perm[i/3]*3+1] = image.data[i+1];
        IMGPerm[perm[i/3]*3+2] = image.data[i+2];
    }

    unsigned char *SVbyte = malloc(sizeof(unsigned char) * 3);
    if(SVbyte == NULL)
    {
        printf("Nu se poate aloca memorie pentru array-ul temporar de descompunere in biti!\n");
        exit(EXIT_FAILURE);
    }

    unsigned char *RNDbyte = malloc(sizeof(unsigned char) * 3);
    if(RNDbyte == NULL)
    {
        printf("Nu se poate aloca memorie pentru array-ul temporar de descompunere in biti!\n");
        exit(EXIT_FAILURE);
    }

    //Conversie numar in vector de bytes (fara byte-ul cel mai semnificativ)
    GetBytes(sv, &SVbyte);
    GetBytes(rand[image.size], &RNDbyte);

    //Primul parametru primeste ceilalte 3 parametrii XOR-ati.
    PixelXOR(&(image.data), SVbyte, IMGPerm, RNDbyte);

    unsigned char *aux;
    for(int index = 3; index < image.size*3; index=index+3)
    {
        GetBytes(rand[image.size+(index/3)], &RNDbyte);
        aux = image.data + index;
        PixelXOR(&aux, aux-3, IMGPerm + index, RNDbyte);
    }

    ImgDatatoBMP(out_filename, image);

    free(SVbyte);
    free(RNDbyte);
    free(perm);
    free(rand);
    free(IMGPerm);
    free(image.data);
    free(image.header);
}

void CalculeazaHistograma(char const in_filename[])
{
    ImgData image;
    BMPtoImgData(in_filename, &image);

    int i;
    int **freq = (int**) malloc(sizeof(int*) * 3);
    freq[0] = (int*) malloc(sizeof(int) * 256);
    freq[1] = (int*) malloc(sizeof(int) * 256);
    freq[2] = (int*) malloc(sizeof(int) * 256);

    if(freq[0] == NULL || freq[1] == NULL || freq[2] == NULL || *freq == NULL)
    {
        printf("Nu se poate aloca memorie pentru vectorul de frecvente al culorilor.");
        exit(EXIT_FAILURE);
    }
    double sB = 0, sG = 0, sR = 0;
    double f2 = (double) image.size/256;

    for(i = 0; i < 256; i++)
    {
        freq[0][i] = 0;
        freq[1][i] = 0;
        freq[2][i] = 0;
    }
    for(i = 0; i < image.size*3; i=i+3)
    {
        freq[0][image.data[i]] ++;
        freq[1][image.data[i+1]] ++;
        freq[2][image.data[i+2]] ++;
    }
    for(i = 0; i < 256; i++)
    {
        sB = sB + ( ((freq[0][i]-f2) * (freq[0][i]-f2)) / f2);
        sG = sG + ( ((freq[1][i]-f2) * (freq[1][i]-f2)) / f2);
        sR = sR + ( ((freq[2][i]-f2) * (freq[2][i]-f2)) / f2);
    }
    printf("Chi-squared test on RGB channels for %s\n", in_filename);
    printf("R: %.2f\n", sR);
    printf("G: %.2f\n", sG);
    printf("B: %.2f\n", sB);

    free(freq[0]);
    free(freq[1]);
    free(freq[2]);
    free(freq);
    free(image.data);
    free(image.header);
}

void GrayScale(ImgData *img)
{
    int i, j, index;
    unsigned char aux;

    for(i = 0; i < img->H; i++)
    {
        for(j = 0; j < img->W*3; j = j+3)
        {
            index = i*(img->W*3)+j;
            aux = 0.299*img->data[index+2] + 0.587*img->data[index+1] + 0.114*img->data[index];
            img->data[index] = img->data[index+1] = img->data[index+2] = aux;
        }
    }
}

double MedieIntensitatiSablon(ImgData img)
{
    int i, j, index;
    unsigned char S;
    double SM, Saux;
    Saux = 0;
    for(i = 0; i < img.H; i++)
    {
        for(j = 0; j < img.W; j++)
        {
            index = i*img.W*3 + j*3;
            S = img.data[index]; //intensitatea pixelului (i,j) al sablonului
            Saux = Saux + S;//suma intensitatilor
        }
    }

    SM = Saux/img.size; //Media intensitatilor
    return SM;
}

double DeviatieIntensitatiSablon(ImgData img)
{
    int i, j, index;
    unsigned char S;
    double SM, dev, Saux;
    Saux = 0;
    SM = MedieIntensitatiSablon(img);
    for(i = 0; i < img.H; i++)
    {
        for(j = 0; j < img.W; j++)
        {
            index = i*img.W*3 + j*3;
            S = img.data[index]; //intensitatea pixelului (i,j) al sablonului
            Saux = Saux + ( (S - SM) * (S - SM) ); //la patrat aici sau dupa for??
        }
    }
    Saux = Saux /(img.size-1);
    dev = sqrt(Saux); //Deviatia intensitatilor pixelilor
    return dev;
}

double MedieIntensitatiImagine(ImgData image, ImgData tmp, int n, int m)
{
    int i, j, index;
    unsigned char fI;
    double fIM, fIaux;
    fIaux = 0;

    for(i = 0; i < tmp.H; i++)
    {
        for(j = 0; j < tmp.W; j++)
        {
            index = (n+i)*image.W*3 + (m+j)*3; //(linie imagine + linie tmp) * latime imagine + coloana imagine + coloana tmp
            if(n+i >= image.H || m+j >= image.W)//a iesit din imagine
            {
                fI = 0;
            }
            else
                fI = image.data[index];
            fIaux = fIaux + fI; // suma intensitatilor
        }
    }
    fIM = fIaux/tmp.size;//Media intensitatilor;
    return fIM;
}

double DeviatieIntensitatiImagine(ImgData image, ImgData tmp, int n, int m)
{
    int i, j, index;
    unsigned char fI;
    double fIM, dev, fIaux;
    fIM = MedieIntensitatiImagine(image, tmp, n, m);

    for(i = 0; i < tmp.H; i++)
    {
        for(j = 0; j < tmp.W; j++)
        {
            index = (n+i)*image.W*3 + (m+j)*3; //(linie imagine + linie tmp) * latime imagine + coloana imagine + coloana tmp
            if(n+i >= image.H || m+j >= image.W)//a iesit din imagine
            {
                fI = 0;
            }
            else
                fI = image.data[index];

            fIaux = fIaux + ( (fI - fIM) * (fI - fIM) ); // suma intensitatilor
        }
    }
    fIaux = fIaux / (tmp.size-1);
    dev = sqrt(fIaux);
    return dev;
}

void draw(ImgData *image, ImgData tmp, int x, int y, char R, char G, char B)
{
    int i, index;
    for(i = 0; i < tmp.W; i++) // linia de sus
    {
        index = y*image->W*3 + (x+i)*3;
        if(index < image->size*3)
        {
            image->data[index] = B;
            image->data[index+1] = G;
            image->data[index+2] = R;
        }

        index = (y+tmp.H-1)*image->W*3 + (x+i)*3;
        if(index < image->size*3)
        {
            image->data[index] = B;
            image->data[index+1] = G;
            image->data[index+2] = R;
        }
    }
    for(i = 0; i < tmp.H; i++)
    {
        index = (y+i)*image->W*3 + (x+tmp.W-1)*3;
        if(index < image->size*3)
        {
            image->data[index] = B;
            image->data[index+1] = G;
            image->data[index+2] = R;
        }

        index = (y+i)*image->W*3 + x*3;
        if(index < image->size*3)
        {
            image->data[index] = B;
            image->data[index+1] = G;
            image->data[index+2] = R;
        }
    }
}

void TemplateMatch(matchData **matches, int *matchesCount, int cifra, char const imagine_filename[], char const tmp_filename[], double prag)
{
    int i, j, x, k, y, index;
    unsigned char S, fI;
    double SM, fIM, devTMP, devIMG, corr, sum;
    int c1;

    ImgData image, tmp;
    BMPtoImgData(imagine_filename, &image);
    BMPtoImgData(tmp_filename, &tmp);
    GrayScale(&image);
    GrayScale(&tmp);

//sablon

    SM = MedieIntensitatiSablon(tmp);
    devTMP = DeviatieIntensitatiSablon(tmp);

//Imagine
    k = 0;
    for(y = 0; y < image.H; y++)//inaltime
    {
        for(x = 0; x < image.W; x++)//latime
        {
            fIM = MedieIntensitatiImagine(image, tmp, y, x);
            devIMG = DeviatieIntensitatiImagine(image, tmp, y, x);
            c1 = (devIMG * devTMP);
            sum = 0; //sum = suma din forumula corelatiei pentru pozitia (x, y) din imagine
            for(i = 0; i < tmp.H; i++)
            {
                for(j = 0; j < tmp.W; j++)
                {
                    index = (y+i)*image.W*3 + (x+j)*3; //(linie imagine + linie tmp) * latime imagine + coloana imagine + coloana tmp
                    if(y+i >= image.H || x+j >= image.W)//a iesit din imagine
                    {
                        fI = 0;
                    }
                    else
                        fI = image.data[index]; //fI intensitatea pixelului din imginea I peste care se aplica sablonul

                    S = tmp.data[i*tmp.W*3 + j*3];//Intensitatea pixelului din sablon

                    sum = sum + (( fI - fIM ) * ( S - SM ))/c1;
                }
            }
            corr = sum/tmp.size;
            if(corr > prag)
            {
                (*matches)[k].cifra = cifra;
                (*matches)[k].x = x;
                (*matches)[k].y = y;
                (*matches)[k].corr = corr;
                k++;
            }
        }
    }
    *matchesCount = k;
    printf("Cifra %d a fost procesata!\n", cifra);
    free(image.data);
    free(image.header);
    free(tmp.data);
    free(tmp.header);
}

int cmp(void const *A, void const *B)
{
    double a = ((matchData*)A)->corr;
    double b = ((matchData*)B)->corr;
    if(a > b)
        return -1;
    if(a < b)
        return 1;
    return 0;
}

int min(int a, int b)
{
    if(a <= b)
        return a;
    return b;
}

int max(int a, int b)
{
    if(a >= b)
        return a;
    return b;
}

int intersectie(matchData m, matchData n, ImgData tmp)
{
    int latime, inaltime;

    //Avem nevoie doar de colturile stanga-sus si dreapta-jos
    dreptunghi md, nd;
    md.stanga = m.x;
    md.sus = m.y;
    md.dreapta = md.stanga + tmp.W;
    md.jos = md.sus + tmp.H;

    nd.stanga = n.x;
    nd.sus = n.y;
    nd.dreapta = nd.stanga + tmp.W;
    nd.jos = nd.sus + tmp.H;

    latime = max(0, (min(md.dreapta, nd.dreapta) - max(md.stanga, nd.stanga)));
    inaltime = max(0, (min(md.jos, nd.jos) - max(md.sus, nd.sus)));
    if((min(md.dreapta, nd.dreapta) - max(md.stanga, nd.stanga)) == 0) latime = 1;
    if((min(md.jos, nd.jos) - max(md.sus, nd.sus)) == 0) inaltime = 1;

    return latime*inaltime;
}

double suprapunere(matchData m1, matchData m2, ImgData tmp)
{
    int inters = intersectie(m1, m2, tmp);
    double x = (double)inters / (double)(tmp.size * 2 - inters);
    return x;
}

void SortDesc(matchData **matches, int count)
{
    qsort(*matches, count, sizeof(matchData), cmp);
}

void RemoveNonMax(matchData **matches, int count, ImgData tmp, double prag)
{
    SortDesc(matches, count);
    int i, j;
    for(i = 0; i < count-1; i++)
    {
        for(j = i+1; j < count; j++)
        {
            if(suprapunere((*matches)[i], (*matches)[j], tmp) > prag && (*matches)[i].x != -1 && (*matches)[j].x != -1)
            {
                (*matches)[j].x = -1; //nu le elimina din vector dar le marcheaza
            }
        }
    }
}

int main()
{
//Criptare/Decriptare

    char *imagineColor, *imagineCriptata, *imagineDecriptata, *cheieSecreta, *sablonCifra;
    imagineColor = (char*) malloc(sizeof(char)*LUNGIME_NUME_FISIER);
    imagineCriptata = (char*) malloc(sizeof(char)*LUNGIME_NUME_FISIER);
    imagineDecriptata = (char*) malloc(sizeof(char)*LUNGIME_NUME_FISIER);
    cheieSecreta = (char*) malloc(sizeof(char)*LUNGIME_NUME_FISIER);
    sablonCifra = (char*) malloc(sizeof(char)*LUNGIME_NUME_FISIER);
    char *imagineTemplate = (char*) malloc(sizeof(char)*LUNGIME_NUME_FISIER);

    if(imagineColor == NULL || imagineCriptata == NULL || imagineDecriptata == NULL || cheieSecreta == NULL || sablonCifra == NULL)
    {
        printf("Nu se poate aloca memorie pentru caile imaginilor!\n");
        exit(EXIT_FAILURE);
    }

    FILE *f = fopen("cripteaza.txt", "r");
    if(f == NULL)
    {
        printf("Eroare la deschiderea fisierului imagini.txt!\n");
        exit(EXIT_FAILURE);
    }
    fscanf(f, "%s %s %s", imagineColor, imagineCriptata, cheieSecreta);
    Cripteaza(imagineColor, imagineCriptata, cheieSecreta);

    f = fopen("decripteaza.txt", "r");
    if(f == NULL)
    {
        printf("Eroare la deschiderea fisierului imagini.txt!\n");
        exit(EXIT_FAILURE);
    }
    fscanf(f, "%s %s %s", imagineCriptata, imagineDecriptata, cheieSecreta);
    Decripteaza(imagineCriptata, imagineDecriptata, cheieSecreta);

    CalculeazaHistograma(imagineColor);
    CalculeazaHistograma(imagineCriptata);

//Template Matching

    int i, nr_sabloane;

    f = fopen("templatematching.txt", "r");
    if(f == NULL)
    {
        printf("Eroare la deschiderea fisierului imagini.txt!\n");
        exit(EXIT_FAILURE);
    }
    fscanf(f, "%s %d ", imagineColor, &nr_sabloane);

    ImgData image, tmp;
    BMPtoImgData(imagineColor, &image);

    matchData *matches;
    matches = (matchData*) malloc(sizeof(matchData) * nr_sabloane * image.size);//alocare spatiu maxim
    if(matches == NULL)
    {
        printf("Nu se poate aloca memorie pentru array-ul de detectii!\n");
        exit(EXIT_FAILURE);
    }

    int *matchesCount;
    int totalCount = 0;
    matchesCount = (int*) malloc(sizeof(int)*nr_sabloane); //10 cifre
    if(matches == NULL)
    {
        printf("Nu se poate aloca memorie pentru array-ul de lungimi asociat array-ului de detectii!\n");
        exit(EXIT_FAILURE);
    }

    for(i = 0; i < nr_sabloane; i++)
    {
        fscanf(f, "%s ", sablonCifra);
        matchData *aux = matches + totalCount;
        TemplateMatch(&aux, &(matchesCount[i]), i, imagineColor, sablonCifra, 0.50);
        totalCount = totalCount + matchesCount[i];
    }

    fscanf(f, "%s", imagineTemplate);

    realloc(matchesCount, totalCount*sizeof(matchData));
    BMPtoImgData(sablonCifra, &tmp);

    const unsigned char colors[10][3] =
    {
        {255, 0, 0}, //rosu
        {255, 255, 0}, //galben
        {0, 255, 0}, //verde
        {0, 255, 255}, //cyan
        {255, 0, 255}, //magenta
        {0, 0, 255}, //albastru
        {192,192, 192}, //argintiu
        {255, 140, 0}, //albastru
        {128, 0, 128}, //magenta
        {128, 0, 0} //albastru
    };

    RemoveNonMax(&matches, totalCount, tmp, 0.20);

    for(i = 0; i < totalCount; i++)
    {
        if(matches[i].x != -1) //e non-maxim
        {
            draw(&image, tmp, matches[i].x, matches[i].y, colors[matches[i].cifra][0], colors[matches[i].cifra][1], colors[matches[i].cifra][2]);
        }
    }

    ImgDatatoBMP(imagineTemplate, image);
    free(matches);
    free(matchesCount);
    free(imagineColor);
    free(imagineCriptata);
    free(imagineDecriptata);
    free(sablonCifra);
    free(cheieSecreta);
    free(image.data);
    free(image.header);
    free(tmp.data);
    free(tmp.header);
}
