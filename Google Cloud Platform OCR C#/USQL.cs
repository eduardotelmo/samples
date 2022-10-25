///////////////////////////////////////////////////////////////////////////////////////////////////////////
// USQL: UnStructured Query Language
///////////////////////////////////////////////////////////////////////////////////////////////////////////

using System;
using System.Text;
using System.Linq;
using System.Collections.Generic;
using System.Globalization;
using System.Text.RegularExpressions;

namespace Google_OCR
{
	public class USQL
	{
		///***************************************************
		/// USQL: UnStructured Query Language
		/// Eduardo Telmo Fonseca Santos
		/// 2017
		///***************************************************

		int min_match;

		public void set_min_match(int _min_match)
		{
			min_match = _min_match;
		}

		public string execute_USQL(string input, string txt)
		{
            string pattern = @"(?i)extrair\s(?'tipo'[a-z]+)(\sentre\s\""(?'termo1'.+)\""\se\s\""(?'termo2'.+)\""((\s)+(maxpalavras(\s)*=(\s)*(?'maxwords'\d+)))?((\s)*(de(\s)+(?'numword1'\d+)(\s)+a(\s)+(?'numword2'\d+)))?)?(\n)?";

			string output = "";

            foreach (Match m in Regex.Matches(input, pattern))
            {
                string termo1 = m.Groups["termo1"].Value;
                string termo2 = m.Groups["termo2"].Value;
                string tipo = m.Groups["tipo"].Value;

                string strmaxwords = m.Groups["maxwords"].Value;
                int maxwords;
                if (strmaxwords != "")
                    maxwords = Int32.Parse(strmaxwords);
                else
                    maxwords = 0;

                string strnumword1 = m.Groups["numword1"].Value;
                int numword1;
                if (strnumword1 != "")
                    numword1 = Int32.Parse(strnumword1);
                else
                    numword1 = 0;

                string strnumword2 = m.Groups["numword2"].Value;
                int numword2;
                if (strnumword2 != "")
                    numword2 = Int32.Parse(strnumword2);
                else
                    numword2 = 0;

                // txtBusca1.Text = termo1;
                // txtBusca2.Text = termo2;


                // Restringe tipo buscado com expressoes regulares
                tipo = tipo.ToLower();
                string type_pattern;
                string match_result = "";
                string nextwords = "";

                // Busca entre os termos ou busca recorte do texto
                if (tipo == "categoria")
                {
                    int txt_length = txt.Length;
                    int txt_crop = 500; // n primeiros caracteres

                    if (txt_crop >= txt_length)
                        txt_crop = txt_length-1;

                    nextwords = txt.Substring(0, txt_crop); // Analisa recorte do texto

					RegexOptions options = RegexOptions.IgnoreCase;

                    type_pattern = @"(escritura)|(contrato)";
					match_result = Regex.Match(nextwords, type_pattern, options).ToString();

					if (match_result == "")
						match_result = "<MISMATCH>";
                }
                else
                {
                    nextwords = search_between_sentences(termo1, termo2, txt, maxwords);
                }

				if (tipo == "palavras")
					match_result = nextwords; // Default

				if (tipo == "data")
				{
					type_pattern = @"[0-9]+ / [0-9]+ / [0-9]+";
					//type_pattern = @"(0[1-9]|[12][0-9]|3[01]) / (0[1-9]|1[0-2]) / ([0-9]*[0-9]{2})";
					match_result = Regex.Match(nextwords, type_pattern).ToString();

					// Remove espacos
					match_result = Regex.Replace(match_result, @"\s+", "", RegexOptions.Multiline);

					if (match_result == "")
						match_result = "<MISMATCH>";
				}

				if (tipo == "maiusculas")
				{
					type_pattern = @"(\b[^a-z\s]+[^a-z]+)+";
					match_result = Regex.Match(nextwords, type_pattern).ToString();

					// Remove pontuacao de separacao entre palavras
					// match_result = Regex.Replace(match_result, @"(\s[.,:;])", "", RegexOptions.Multiline);

					if (match_result == "")
						match_result = "<MISMATCH maiusculas>";
				}

				if (tipo == "minusculas")
				{
					type_pattern = @"(\b[^A-Z\s]+[^A-Z]+)+";
					match_result = Regex.Match(nextwords, type_pattern).ToString();

					// Remove pontuacao de separacao entre palavras
					// match_result = Regex.Replace(match_result, @"(\s[.,:;])", "", RegexOptions.Multiline);

					if (match_result == "")
						match_result = "<MISMATCH minusculas>";
				}

				if (tipo == "cpf")
				{
					type_pattern = @"[0-9]{3}[\s]?\.[\s]?[0-9]{3}[\s]?\.[\s]?[0-9]{3}[\s]?-[\s]?[0-9]{2}";
					match_result = Regex.Match(nextwords, type_pattern).ToString();

					// Remove espacos
					match_result = Regex.Replace(match_result, @"\s+", "", RegexOptions.Multiline);

					if (match_result == "")
						match_result = "<MISMATCH>";
				}

				if (tipo == "rg")
				{
					type_pattern = @"([0-9]{2}[\s]?\.[\s]?)?[0-9]{3}[\s]?\.[\s]?[0-9]{3}[\s]?-[\s]?[0-9]{2}";
					match_result = Regex.Match(nextwords, type_pattern).ToString();

					// Remove espacos
					match_result = Regex.Replace(match_result, @"\s+", "", RegexOptions.Multiline);

					if (match_result == "")
						match_result = "<MISMATCH>";
				}

				if (tipo == "expedidor")
				{
					type_pattern = @"[A-Z]+ / [A-Z]{2}";
					match_result = Regex.Match(nextwords, type_pattern).ToString();

					// Remove espacos
					match_result = Regex.Replace(match_result, @"\s+", "", RegexOptions.Multiline);

					if (match_result == "")
						match_result = "<MISMATCH>";
				}

				if (tipo == "pfpj")
				{
					// Verifica se tem RG
					type_pattern = @"([0-9]{2}[\s]?\.[\s]?)?[0-9]{3}[\s]?\.[\s]?[0-9]{3}[\s]?-[\s]?[0-9]{2}";
					match_result = Regex.Match(nextwords, type_pattern).ToString();

					// Remove espacos
					match_result = Regex.Replace(match_result, @"\s+", "", RegexOptions.Multiline);

					if (match_result == "")
						match_result = "PJ";
					else
						match_result = "PF";
				}

				if (tipo == "estadocivil")
				{
					RegexOptions options = RegexOptions.IgnoreCase;

					type_pattern = @"(casad[ao])|(solteir[ao])|(divorciad[ao])|(vi.v[ao])|(uni.o\sest.vel)";
					match_result = Regex.Match(nextwords, type_pattern, options).ToString();

					if (match_result == "")
						match_result = "<MISMATCH>";
				}

				if (tipo == "valor")
				{
					RegexOptions options = RegexOptions.IgnoreCase;

					type_pattern = @"[A-Z]*\$ .+ , [0-9]{2}";
					match_result = Regex.Match(nextwords, type_pattern, options).ToString();

					// Remove espacos
					match_result = Regex.Replace(match_result, @"\s+", "", RegexOptions.Multiline);

					if (match_result == "")
						match_result = "<MISMATCH>";
				}

                // Extrai palavras pelo indice
				if (numword1 != 0)
				{
                    int list_length;

					List<string> list_words = match_result.Split(' ').ToList();
					match_result = "";

					// Verificacao de limites dos indices
					list_length = list_words.Count;
					if (numword1 > list_length)
						numword1 = list_length;
					if (numword2 > list_length)
						numword2 = list_length;
					for (int i = numword1 - 1; i < numword2; i++)
						match_result += list_words[i] + " ";
				}

				// !!! txtQueryResults.Text += match_result + "\r\n";

				output += match_result + "\r\n";
			}

			return output;

		}

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////


		public void search_after_sentence(string str1, string txt)
		{
			string found = "";
			int matchindex;

			string nextwords = ClosestMatch(str1, txt, out found, out matchindex, min_match); // Busca str1 em txt

			// !!! form.txtEncontrou.Text = found;
			// !!! txtProximas.Text = nextwords;
		}

		///////////////////////////////////////////////////////////////////////////////////////////////////////////


		public string search_between_sentences(string str1, string str2, string txt, int maxwords = 0)
		{
			string found = "";
			int matchindex;

			string nextwords = ClosestMatch(str1, txt, out found, out matchindex, min_match); // Busca str1 em txt

			// !!! txtEncontrou.Text = found;
			// !!! txtProximas.Text = nextwords;

			//////////////////////////////////////////////////////////////////////////////

			if (found != "")
			{
				nextwords = ClosestMatch(str2, txt, out found, out matchindex, min_match, matchindex + 1, maxwords); // Busca str2 em txt

				// !!! txtEncontrou.Text = found;
				// !!! txtProximas.Text = nextwords;
			}

			return nextwords;
		}

		///////////////////////////////////////////////////////////////////////////////////////////////////////////

		public string ClosestMatch(string term, string terms, out string found, out int matchindex, int min_match = 0, int startindex = 0, int maxwords = 0) // !!!
		{
			Levenshtein lev = new Levenshtein();
			int maxj;

			// Cria uma lista com as palavras separadas por espacos
			List<string> list_term = term.Split(' ').ToList();
			List<string> list_terms = terms.Split(' ').ToList();

			// Inicializa vetor com valor maximo de inteiro
			int[] counter = Enumerable.Repeat(Int32.MaxValue, list_terms.Count).ToArray();

			// Calcula indice maximo a partir do limite para o numero maximo de palavras
			int endindex;
			if (maxwords != 0)
				endindex = startindex + maxwords + 2;
			else
				endindex = list_terms.Count;

			int num_words_term = list_term.Count;

			// Loop de comparacao de termos no texto
			for (int i = startindex; i < endindex; i++)
			{
				string str1;
				string str2;

				// Comparacao de multiplas palavras
				string terms_cropped = list_terms[i];
				if (num_words_term > 1)
				{
					maxj = list_terms.Count;
					for (int j = i + 1; (j < i + num_words_term) && (j < maxj); j++)
						terms_cropped = terms_cropped + " " + list_terms[j];
				}
				str1 = NormalizeWord(term);
				str2 = NormalizeWord(terms_cropped);

				// Calcula distancia de Levenshtein entre as palavras (distancia minima de edicao)
				counter[i] = lev.iLD(str1, str2); // Mais eficiente
			}

			// Encontra termos com a menor distancia de Levenshtein
			int min = counter.Min();
			int index = counter.ToList().FindIndex(t => t == min);
			matchindex = index + num_words_term - 1;
			found = list_terms[index];

			string ret = "";

			if (startindex > 0)
			{
				if ((maxwords == 0) || ((matchindex - startindex) <= maxwords))
				{
					for (int i = startindex; i < matchindex - num_words_term + 1; i++)
						ret += list_terms[i] + " ";
				}
				else
					ret = "<MAXWORDS>"; // !!! Excedeu o numero maximo de palavaras
			}
			else
			{
				// Exibe as proximas palavras depois dos termos encontrados
				maxj = list_terms.Count;
				int dj = maxwords;
				if (maxwords == 0) // Para testes
					dj = 3;
				for (int j = 0; j < dj; j++)
				{
					int k = index + j + list_term.Count;
					if (k < maxj)
						ret += list_terms[k] + " ";
					else
						break;
				}
			}

			// Calcula percentual de semelhanca com os termos buscados
			float perc_match;
			perc_match = 100.0f * (1.0f - (float)min / (float)term.Length);
			if (perc_match < min_match)
			{
				ret = "<MISMATCH>";
				found = "";
			}
			return ret;
		}

		///////////////////////////////////////////////////////////////////////////////////////////////////////////

		public string NormalizeText(string text)
		{
			// Separa simbolos das palavras e numeros, colocando espacos antes e depois
			text = Regex.Replace(text, @"([\.,/:;()-])", @" $1 ", RegexOptions.Multiline);

			// Substitui multiplos espacos por espaco simples
			text = Regex.Replace(text, @"\s+", " ", RegexOptions.Multiline);
			// text = Regex.Replace(text, @"\s+", " ");

			return text;
		}

		///////////////////////////////////////////////////////////////////////////////////////////////////////////

		public string NormalizeWord(string text)
		{
			// Modifica entradas para lowercase sem acentos
			text = text.ToLower();

			// Remove simbolos (pode causar problemas em alguns casos)
			//!!! text = Regex.Replace(text, "[.,;:]", "");

			// Remove acentos
			StringBuilder sbReturn = new StringBuilder();
			var arrayText = text.Normalize(NormalizationForm.FormD).ToCharArray();
			foreach (char letter in arrayText)
			{
				if (CharUnicodeInfo.GetUnicodeCategory(letter) != UnicodeCategory.NonSpacingMark)
					sbReturn.Append(letter);
			}
			return sbReturn.ToString();
		}


	}

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
// Distancia de Levenshtein (codigo externo)
///////////////////////////////////////////////////////////////////////////////////////////////////////////
public class Levenshtein
{
	///***************************************************
	/// Calcula distancia de Levenshtein entre strings
	/// Versao eficiente em memoria
	/// Metodo sobrecarregado (absoluto e percentual)
	///***************************************************

	public int iLD(String sRow, String sCol)
	{
		int perc_diff;
		return iLD_core(sRow, sCol, out perc_diff);
	}

	public int iLD(String sRow, String sCol, out int perc_diff)
	{
		return iLD_core(sRow, sCol, out perc_diff);
	}

	public int iLD_core(String sRow, String sCol, out int perc_diff)
	{
		int RowLen = sRow.Length;  // comprimento de sRow
		int ColLen = sCol.Length;  // comprimento de sCol
		int RowIdx;                // indice de iteracao atraves de sRow
		int ColIdx;                // indice de iteracao atraves de sCol
		char Row_i;                // i-esimo caracter de sRow
		char Col_j;                // j-esimo caracter de sCol
		int cost;                   // custo

		perc_diff = 100; // Caso saia antes da comparacao

		/// Testa comprimento da string
		if (Math.Max(sRow.Length, sCol.Length) > Math.Pow(2, 31))
			throw (new Exception("\nMaximo comprimento da string em Levenshtein.iLD é " + Math.Pow(2, 31) + ".\nO seu é " + Math.Max(sRow.Length, sCol.Length) + "."));

		// Passo 1

		if (RowLen == 0)
		{
			return ColLen;
		}

		if (ColLen == 0)
		{
			return RowLen;
		}

		/// Cria os dois vetores
		int[] v0 = new int[RowLen + 1];
		int[] v1 = new int[RowLen + 1];
		int[] vTmp;



		/// Passo 2
		/// Inicializa o primeiro vetor
		for (RowIdx = 1; RowIdx <= RowLen; RowIdx++)
		{
			v0[RowIdx] = RowIdx;
		}

		// Passo 3

		/// Para cada coluna
		for (ColIdx = 1; ColIdx <= ColLen; ColIdx++)
		{
			/// Ajusta o elemento 0 para o numero da coluna
			v1[0] = ColIdx;

			Col_j = sCol[ColIdx - 1];


			// Passo 4

			/// Para cada linha
			for (RowIdx = 1; RowIdx <= RowLen; RowIdx++)
			{
				Row_i = sRow[RowIdx - 1];


				// Passo 5

				if (Row_i == Col_j)
				{
					cost = 0;
				}
				else
				{
					cost = 1;
				}

				// Passo 6

				/// Busca o minimo
				int m_min = v0[RowIdx] + 1;
				int b = v1[RowIdx - 1] + 1;
				int c = v0[RowIdx - 1] + cost;

				if (b < m_min)
				{
					m_min = b;
				}
				if (c < m_min)
				{
					m_min = c;
				}

				v1[RowIdx] = m_min;
			}

			/// Troca os vetores
			vTmp = v0;
			v0 = v1;
			v1 = vTmp;

		}


		// Passo 7

		/// Valor entre 0 - 100
		/// 0==casamento perfeito 100==totalmente diferente
		/// 
		/// Os vetores sao trocados uma ultima vez no final do ultimo loop,
		/// por isto o resultado esta agora em v0 ao inves de v1
		//System.Console.WriteLine("Levenshein iDist=" + v0[RowLen]);
		int max = System.Math.Max(RowLen, ColLen);
		perc_diff = ((100 * v0[RowLen]) / max); // Saida percentual (normalizada)

		return v0[RowLen]; // Retorna o numero de operacoes para transformar uma string na outra (diferenca)
	}

	///*****************************
	/// Calcula o minimo
	///*****************************

	private int Minimum(int a, int b, int c)
	{
		int mi = a;

		if (b < mi)
		{
			mi = b;
		}
		if (c < mi)
		{
			mi = c;
		}

		return mi;
	}

	///*************************************************
	/// Calcula a distancia de Levenshtein  
	/// Versao sem otimizacao de memoria    
	/// Metodo sobrecarregado (absoluto e percentual)   
	///*************************************************

	public int LD(String sRow, String sCol)
	{
		int perc_diff;
		return LD_core(sRow, sCol, out perc_diff);
	}

	public int LD(String sRow, String sCol, out int perc_diff)
	{
		return LD_core(sRow, sCol, out perc_diff);
	}

	public int LD_core(String sNew, String sOld, out int perc_diff)
	{
		int[,] matrix;              // matriz
		int sNewLen = sNew.Length;  // comprimento de sNew
		int sOldLen = sOld.Length;  // comprimento de sOld
		int sNewIdx; // indice para iteracao atraves de sNew
		int sOldIdx; // indice para iteracao atraves de sOld
		char sNew_i; // i-esimo caracter de sNew
		char sOld_j; // j-esimo caracter de sOld
		int cost; // custo

		perc_diff = 100; // Caso saia antes da comparacao

		/// Testa comprimento da string
		if (Math.Max(sNew.Length, sOld.Length) > Math.Pow(2, 31))
			throw (new Exception("\nO tamanho maximo da string em Levenshtein.LD é " + Math.Pow(2, 31) + ".\nO seu é " + Math.Max(sNew.Length, sOld.Length) + "."));

		// Passo 1

		if (sNewLen == 0)
		{
			return sOldLen;
		}

		if (sOldLen == 0)
		{
			return sNewLen;
		}

		matrix = new int[sNewLen + 1, sOldLen + 1];

		// Passo 2

		for (sNewIdx = 0; sNewIdx <= sNewLen; sNewIdx++)
		{
			matrix[sNewIdx, 0] = sNewIdx;
		}

		for (sOldIdx = 0; sOldIdx <= sOldLen; sOldIdx++)
		{
			matrix[0, sOldIdx] = sOldIdx;
		}

		// Passo 3

		for (sNewIdx = 1; sNewIdx <= sNewLen; sNewIdx++)
		{
			sNew_i = sNew[sNewIdx - 1];

			// Passo 4

			for (sOldIdx = 1; sOldIdx <= sOldLen; sOldIdx++)
			{
				sOld_j = sOld[sOldIdx - 1];

				// Passo 5

				if (sNew_i == sOld_j)
				{
					cost = 0;
				}
				else
				{
					cost = 1;
				}

				// Passo 6

				matrix[sNewIdx, sOldIdx] = Minimum(matrix[sNewIdx - 1, sOldIdx] + 1, matrix[sNewIdx, sOldIdx - 1] + 1, matrix[sNewIdx - 1, sOldIdx - 1] + cost);

			}
		}

		// Passo 7

		/// Valor entre 0 - 100
		/// 0==casamento perfeito 100==totalmente diferente
		//System.Console.WriteLine("Dist=" + matrix[sNewLen, sOldLen]);
		int max = System.Math.Max(sNewLen, sOldLen);
		perc_diff = (100 * matrix[sNewLen, sOldLen]) / max; // Saida percentual (normalizada)
		return matrix[sNewLen, sOldLen]; // Retorna o numero de operacoes para transformar uma string na outra (diferenca)
	}
}
