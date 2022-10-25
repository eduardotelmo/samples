///////////////////////////////////////////////////////////////////////////////////////////////////////////
// Signature_Detector: Detecta assinatura em imagem utilizando OpenCV
///////////////////////////////////////////////////////////////////////////////////////////////////////////

using Emgu.CV;
using Emgu.CV.CvEnum;
using Emgu.CV.Structure;

using System;
using System.Drawing;

namespace Google_OCR
{
    public class Signature_Detector
    {
        // Containers das imagens para deteccao de assinatura
        private static Mat c = new Mat(); // Imagem reduzida
        private static Mat d = new Mat(); // Imagem sharpen
        private static Mat dd = new Mat(); // Retangulo contendo apenas a maior imagem conectada
        private static Mat p1 = new Mat(); // Imagem Prewitt vertical 1
        private static Mat p2 = new Mat(); // Imagem Prewitt vertical 2
        private static Mat cb = new Mat(); // Bloco para analise de assinatura

        // Variaveis do algoritmo para detectar a maior curva de um recorte da imagem
        // Traduzido do Matlab com adaptacoes para C#
        private static int[] iib;
        private static int[] jjb;
        private static int[] best_iib;
        private static int[] best_jjb;
        private static int npb;
        private static int best_npb;
        private static int mmb;
        private static int nnb;
        private static int wi;
        private static int wj;

        public int Detect_Signature(string imgFile)
        {
            // Para TESTES
            // Page16.jpg
            // string imgFile = @"C:\Users\Esmae\Desktop\Page16.jpg";
            // string imgFile = @"\\Mac\Home\Desktop\Page16.jpg";

            Mat picture = new Mat(imgFile); // Le arquivo de imagem
            // CvInvoke.Imshow("Imagem original", picture); // Abre janela com imagem

            // Reduz a imagem 4x
            int w, h;
            w = picture.Width / 4;
            h = picture.Height / 4;
            CvInvoke.Resize(picture, c, new Size(w, h), 0, 0, Inter.Linear);
            // CvInvoke.Imshow("Imagem reduzida 4x", c); // Abre janela com imagem

            // Converte em grayscale
            CvInvoke.CvtColor(c, c, ColorConversion.Bgr2Gray);

            // Filtro Sharpen
            Matrix<float> sharpen_kernel = new Matrix<float>(new float[3, 3] { { 0, -1, 0 }, { -1, 5, -1 }, { 0, -1, 0 } });
            d = new Mat(c.Rows, c.Cols, DepthType.Cv8U, 1);
            CvInvoke.Filter2D(c, d, sharpen_kernel, new Point(-1, -1));
            // CvInvoke.Imshow("Imagem Sharpen", d); // Abre janela com imagem

            // Remove vertical edges (Prewitt 1)
            p1 = new Mat(d.Rows, d.Cols, DepthType.Cv8U, 1);
            d.CopyTo(p1);
            // Negativo da imagem
            CvInvoke.BitwiseNot(p1, p1);
            // Filtro Prewitt
            Matrix<float> prewitt_kernel1 = new Matrix<float>(new float[7, 3] { { 1, 0, -1 }, { 1, 0, -1 }, { 1, 0, -1 }, { 1, 0, -1 }, { 1, 0, -1 }, { 1, 0, -1 }, { 1, 0, -1 } });
            CvInvoke.Filter2D(p1, p1, prewitt_kernel1, new Point(-1, -1));
            // Negativo da imagem
            CvInvoke.BitwiseNot(p1, p1);
            // CvInvoke.Imshow("Imagem Prewitt Vertical 1", p1); // Abre janela com imagem

            // Remove vertical edges (Prewitt 2)
            p2 = new Mat(d.Rows, d.Cols, DepthType.Cv8U, 1);
            d.CopyTo(p2);
            // Negativo da imagem
            CvInvoke.BitwiseNot(p2, p2);
            // Filtro Prewitt
            Matrix<float> prewitt_kernel2 = new Matrix<float>(new float[7, 3] { { -1, 0, 1 }, { -1, 0, 1 }, { -1, 0, 1 }, { -1, 0, 1 }, { -1, 0, 1 }, { -1, 0, 1 }, { -1, 0, 1 } });
            CvInvoke.Filter2D(p2, p2, prewitt_kernel2, new Point(-1, -1));
            // Negativo da imagem
            CvInvoke.BitwiseNot(p2, p2);
            // CvInvoke.Imshow("Imagem Prewitt Vertical 2", p2); // Abre janela com imagem

            // Cria Matrix para poder acessar via indices
            Matrix<Byte> p1mat = new Matrix<Byte>(p1.Rows, p1.Cols, p1.DataPointer);
            Matrix<Byte> p2mat = new Matrix<Byte>(p2.Rows, p2.Cols, p2.DataPointer);

            // Remocao de linhas longas horizontais e verticais com moving window (analise de blocos da imagem)
            int mm = d.Rows;
            int nn = d.Cols;
            wi = (int)Math.Round((double)mm / 16); // Poderia fixar este valor e o tamanho total da imagem
            wj = wi;

            // Aloca bloco de trabalho
            cb = new Mat(wi, wj, DepthType.Cv8U, 1);

            // Console.WriteLine("Removendo linhas");
            // Nao pode alterar diretamente variaveis MAT mas pode alterar variaveis MATRIX obtidas a partir do MAT
            Matrix<Byte> dmat = new Matrix<Byte>(mm, nn, d.DataPointer);

            // Pre-processamento da imagem (converte em B&W e remove linhas)
            for (int i = 0; i < (mm - wi + 1); i += wi)
            {
                for (int j = 0; j < (nn - wj + 1); j += wj)
                {
                    // Copia bloco da imagem
                    Matrix<byte> dm = new Matrix<byte>(wi, wj);
                    Matrix<byte> dmorg = new Matrix<byte>(wi, wj);
                    for (int ii = 0; ii < wi; ii++)
                        for (int jj = 0; jj < wj; jj++)
                            dm[ii, jj] = dmat[i + ii, j + jj];

                    // Conversao em B&W
                    int md = wi;
                    int nd = wj;

                    double dmean = (double)dm.Sum / (double)(md * nd);
                    double thr = Math.Max(dmean, (double)5);

                    // Converte imagem em preto e branco
                    for (int ii = 0; ii < md; ii++)
                    {
                        for (int jj = 0; jj < nd; jj++)
                        {
                            if (dm[ii, jj] < thr)
                                dm[ii, jj] = 0;
                            if (dm[ii, jj] >= thr)
                                dm[ii, jj] = 255;
                            // Negativo da imagem
                            dm[ii, jj] = (byte)(255 - dm[ii, jj]);
                        }
                    }

                    // Guarda bloco original
                    dm.CopyTo(dmorg);

                    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                    // Remove linhas horizontais
                    int wl = 10;
                    Mat se1 = CvInvoke.GetStructuringElement(ElementShape.Rectangle, new Size(wl, 1), new Point(-1, -1));
                    CvInvoke.Erode(dm, dm, se1, new Point(-1, -1), 1, BorderType.Constant, CvInvoke.MorphologyDefaultBorderValue);
                    for (int ii = 0; ii < wi; ii++)
                        for (int jj = 0; jj < wj; jj++)
                            dm[ii, jj] = (byte)(dmorg[ii, jj] - dm[ii, jj]);
                    dm.CopyTo(dmorg);

                    // Remove linhas verticais
                    wl = 24;
                    Mat se2 = CvInvoke.GetStructuringElement(ElementShape.Rectangle, new Size(1, wl), new Point(-1, -1));
                    CvInvoke.Erode(dm, dm, se2, new Point(-1, -1), 1, BorderType.Constant, CvInvoke.MorphologyDefaultBorderValue);
                    for (int ii = 0; ii < wi; ii++)
                        for (int jj = 0; jj < wj; jj++)
                            dm[ii, jj] = (byte)(dmorg[ii, jj] - dm[ii, jj]);
                    dm.CopyTo(dmorg);

                    // Positivo da imagem
                    for (int ii = 0; ii < wi; ii++)
                        for (int jj = 0; jj < wj; jj++)
                            dm[ii, jj] = (byte)(255 - dm[ii, jj]);
                    dm.CopyTo(dmorg);

                    // Cola bloco da imagem
                    for (int ii = 0; ii < wi; ii++)
                        for (int jj = 0; jj < wj; jj++)
                            dmat[i + ii, j + jj] = dm[ii, jj];

                }
            }

            // Aloca vetores de coordenadas dos pontos brancos conectados detectados no bloco
            iib = new int[wi * wj];
            jjb = new int[wi * wj];
            best_iib = new int[wi * wj];
            best_jjb = new int[wi * wj];

            // Loop principal de analise dos blocos da imagem
            int cnt = 0;
            for (int i = 0; i < (mm - wi + 1); i += wi)
            {
                for (int j = 0; j < (nn - wj + 1); j += wj)
                {
                    // Copia bloco da imagem (d)
                    Matrix<byte> dm = new Matrix<byte>(wi, wj);
                    int md = dm.Rows;
                    int nd = dm.Cols;
                    for (int ii = 0; ii < wi; ii++)
                        for (int jj = 0; jj < wj; jj++)
                            dm[ii, jj] = dmat[i + ii, j + jj];

                    // Copia bloco da imagem (p1c e p2c)
                    Matrix<byte> p1c = new Matrix<byte>(wi, wj);
                    Matrix<byte> p2c = new Matrix<byte>(wi, wj);
                    for (int ii = 0; ii < wi; ii++)
                        for (int jj = 0; jj < wj; jj++)
                        {
                            p1c[ii, jj] = p1mat[i + ii, j + jj];
                            p2c[ii, jj] = p2mat[i + ii, j + jj];
                        }

                    // Binarizacao para deteccao morfologica da maior componente conectada (verificar se e necessario)
                    for (int ii = 0; ii < wi; ii++)
                        for (int jj = 0; jj < wj; jj++)
                        {
                            if (dm[ii, jj] > 0)
                                dm[ii, jj] = 255;
                            dm[ii, jj] = (byte)(255 - dm[ii, jj]);
                        }

                    // Bloco de trabalho para detectar a maior curva
                    Matrix<byte> cbmat = new Matrix<byte>(wi, wj, cb.DataPointer);
                    dm.CopyTo(cbmat);

                    /////////////////////////////////////////////////////////////////////////////////////////////

                    // Encontra a maior componente conectada
                    mmb = wi;
                    nnb = wj;
                    bwconncomp_main();
                    cb.CopyTo(dm);

                    // Encontra o menor retangulo que contem a componente conectada
                    int minii = 0, minjj = 0, maxii = 0, maxjj = 0;
                    if (best_npb > 0)
                    {
                        minii = best_iib[0];
                        minjj = best_jjb[0];
                        maxii = best_iib[0];
                        maxjj = best_jjb[0];
                    }
                    for (int ii = 0; ii < best_npb; ii++)
                    {
                        if (best_iib[ii] < minii)
                            minii = best_iib[ii];
                        if (best_jjb[ii] < minjj)
                            minjj = best_jjb[ii];
                        if (best_iib[ii] > maxii)
                            maxii = best_iib[ii];
                        if (best_jjb[ii] > maxjj)
                            maxjj = best_jjb[ii];
                    }
                    int wib = Math.Abs(maxii - minii);
                    int wjb = Math.Abs(maxjj - minjj);

                    // Bounding Box
                    // Copia o menor retangulo contendo a maior componente conectada
                    // Imagem com os valores entre 0 e 1
                    dd = new Mat(wib, wjb, DepthType.Cv8U, 1);
                    Matrix<byte> ddmat = new Matrix<byte>(wib, wjb, dd.DataPointer);
                    if ((best_npb > 0) && (wib > 0) && (wjb > 0))
                    {
                        for (int ii = 0; ii < wib; ii++)
                            for (int jj = 0; jj < wjb; jj++)
                            {
                                ddmat[ii, jj] = dm[ii + minii, jj + minjj]; // Copia para dd (a ser analisado)
                                // dmat[i + ii, j + jj] = dm[ii + minii, jj + minjj]; // Copia para a tela

                                // Recorta p1c e p2c
                                p1c[ii, jj] = p1c[ii + minii, jj + minjj]; // Copia para p1c (a ser analisado)
                                p2c[ii, jj] = p1c[ii + minii, jj + minjj]; // Copia para p2c (a ser analisado)

                                // Limita no intervalo 0...1
                                if (ddmat[ii, jj] > 0)
                                    ddmat[ii, jj] = 1;
                                if (p1c[ii, jj] > 0)
                                    p1c[ii, jj] = 1;
                                if (p2c[ii, jj] > 0)
                                    p2c[ii, jj] = 1;
                            }

                        // Para testes de bounding box
                        // CvInvoke.Imshow("Bounding Box " + i*wi+j, ddmat*255); // Open window with image

                        ////////////////////////////////////////////////////////////////////////////////////////

                        // Analise das propriedades estatisticas da Bounding Box (Features da imagem)
                        double dens;
                        dens = 0;
                        for (int ii = 0; ii < wib; ii++)
                            for (int jj = 0; jj < wjb; jj++)
                                dens = dens + ddmat[ii, jj];
                        // dens = dens / (double)(wi * wj);
                        dens = dens / (double)(wib * wjb);
                        // Console.WriteLine(dens);

                        int mid = (int)Math.Round((double)(wib / 2.0));
                        // Console.WriteLine(mid);

                        int border = (int)Math.Round((double)(wib / 8.0));
                        if (border < 1)
                            border = 1;
                        // Console.WriteLine(border);

                        double ud = 0;
                        double ld = 0;
                        double mu, nu;
                        double ml, nl;
                        mu = mid;
                        nu = wjb;
                        ml = Math.Abs(wib - mid);
                        nl = wjb;
                        ud = 0;
                        for (int ii = 0; ii < mid; ii++)
                            for (int jj = 0; jj < nu; jj++)
                                ud = ud + ddmat[ii, jj];
                        ud = ud / (double)(mu * nu);
                        ld = 0;
                        for (int ii = mid; ii < wib; ii++)
                            for (int jj = 0; jj < nu; jj++)
                                ld = ld + ddmat[ii, jj];
                        ld = ld / (double)(ml * nl);

                        double pd = Math.Abs(ud - ld);
                        // Console.WriteLine(pd);

                        int ve, mve;
                        mve = 0;
                        for (int ii = 0; ii < wib; ii++)
                        {
                            ve = 0;
                            for (int jj = 0; jj < wjb; jj++)
                                ve = ve + p1c[ii, jj];
                            if (ve > mve)
                                mve = ve;
                        }
                        // Console.WriteLine(mve);

                        // Identifica se e uma assinatura pelas propriedades estatisticas
                        int handwritten = 0;
                        if ((dens >= 0.07) && (dens <= 0.2) && (mve >= 14) && (pd > 0.01) && (wib > 2) && (wjb > 2))
                        {
                            handwritten = 1;
                            cnt = cnt + 1;
                        }

                        // Original (Matlab):
                        /*
                        if ((dens >= 0.01) && (mve >= 15) && (pd > 0.008) && (wib > 2) && (wjb > 2))
                        {
                            handwritten = 1;
                            cnt = cnt + 1;
                        }
                        */

                        // Imagem com blocos de assinatura identificada exibidos em negativo
                        if (handwritten == 1)
                        {
                            for (int ii = 0; ii < wi; ii++)
                                for (int jj = 0; jj < wj; jj++)
                                    dmat[i + ii, j + jj] = dm[ii, jj];
                        }
                        else
                        {
                            for (int ii = 0; ii < wi; ii++)
                                for (int jj = 0; jj < wj; jj++)
                                    dmat[i + ii, j + jj] = (byte)(255 - dm[ii, jj]);
                        }

                    }
                    else
                    {
                        // Blocos descartados sumariamente
                        for (int ii = 0; ii < wi; ii++)
                            for (int jj = 0; jj < wj; jj++)
                                dmat[i + ii, j + jj] = (byte)(255 - dm[ii, jj]);
                    }

                    ////////////////////////////////////////////////////////////////////////////////////////

                }
            }

            // Desenha grade de retangulos
            for (int i = 0; i < (mm - wi + 1); i += wi)
            {
                for (int j = 0; j < (nn - wj + 1); j += wj)
                {
                    CvInvoke.Rectangle(dmat, new Rectangle(j, i, wj, wi), new MCvScalar(127), 1);
                }
            }


            // Exibe resultado para testes:
            /*
            Console.WriteLine("Blocos com assinatura: " + cnt);

            CvInvoke.Imshow("Imagem processada", d); // Open window with image
            Console.WriteLine("Pressione qualquer tecla sobre uma imagem para sair");
            CvInvoke.WaitKey();
            */

            // Retorna o numero de blocos com assinaturas detectadas
            return cnt;
        }

        //////////////////////////////////////////////////////////////////////////////////////////////////////

        private static void bwconncomp_main()
        {
            // Funcao principal para identificar a maior componente conectada na imagem (linha branca)
            Matrix<byte> cbmat = new Matrix<byte>(wi, wj, cb.DataPointer);

            best_npb = 0;
            for (int i = 0; i < mmb; i++)
                for (int j = 0; j < nnb; j++)
                {
                    npb = 0;
                    if (cbmat[i, j] > 0)
                    {
                        // Chamada recursiva para detectar componentes detectadas a partir do pixel (i,j)
                        bwconncomp(i, j);

                        if (npb > best_npb)
                        {
                            best_npb = npb;
                            // Guarda as coordenadas dos pontos brancos conectados
                            iib.CopyTo(best_iib, 0);
                            jjb.CopyTo(best_jjb, 0);
                        }

                    }

                }

            // Desenha maior componente conectada para testes
            for (int i = 0; i < best_npb; i++)
                cbmat[best_iib[i], best_jjb[i]] = 255;

            // Console.WriteLine(best_npb);
        }

        //////////////////////////////////////////////////////////////////////////////////////////////////////

        private static void bwconncomp(int i, int j)
        {
            // Funcao recursiva para identificar a maior componente conectada a partir do pixel (i,j)
            Matrix<byte> cbmat = new Matrix<byte>(wi, wj, cb.DataPointer);

            cbmat[i, j] = 0;

            // Guarda as coordenadas do ponto branco
            iib[npb] = i;
            jjb[npb] = j;
            npb = npb + 1;

            // Cada pixel tem 8 pixels vizinhos (8 recursoes)

            // Cima / Baixo / Esquerda / Direita
            if ((i - 1) >= 0)
                if (cbmat[i - 1, j] > 0)
                    bwconncomp(i - 1, j);
            if ((j - 1) >= 0)
                if (cbmat[i, j - 1] > 0)
                    bwconncomp(i, j - 1);
            if ((i + 1) < mmb)
                if (cbmat[i + 1, j] > 0)
                    bwconncomp(i + 1, j);
            if ((j + 1) < nnb)
                if (cbmat[i, j + 1] > 0)
                    bwconncomp(i, j + 1);

            // Diagonais
            if ((i - 1) >= 0 && (j - 1) >= 0)
                if (cbmat[i - 1, j - 1] > 0)
                    bwconncomp(i - 1, j - 1);
            if ((i - 1) >= 0 && (j + 1) < nnb)
                if (cbmat[i - 1, j + 1] > 0)
                    bwconncomp(i - 1, j + 1);
            if ((i + 1) < mmb && (j - 1) >= 0)
                if (cbmat[i + 1, j - 1] > 0)
                    bwconncomp(i + 1, j - 1);
            if ((i + 1) < mmb && (j + 1) < nnb)
                if (cbmat[i + 1, j + 1] > 0)
                    bwconncomp(i + 1, j + 1);
        }
    }
}
