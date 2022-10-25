using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using Google.Apis.Auth.OAuth2;
using Google.Apis.Services;
using Google.Apis.Vision.v1;
using Google.Apis.Vision.v1.Data;
using Newtonsoft.Json;

namespace Google_OCR
{
    public class GoogleAnnotate
    {
		public int Page_angle { get; set; }

		public string ApplicationName { get { return "Ocr"; } }
        public string JsonResult { get; set; }
        public string TextResult { get; set; }
        public string Error { get; set; }

		public enum Platform
		{
			Windows,
			Linux,
			Mac
		}

		public string FileSlash = "\\"; // Windows

		public static Platform RunningPlatform()
		{
			switch (Environment.OSVersion.Platform)
			{
				case PlatformID.Unix:
					if (Directory.Exists("/Applications")
						& Directory.Exists("/System")
						& Directory.Exists("/Users")
						& Directory.Exists("/Volumes"))
						return Platform.Mac;
					else
						return Platform.Linux;

				case PlatformID.MacOSX:
					return Platform.Mac;

				default:
					return Platform.Windows;
			}
		}

		private string JsonKeypath
        {
			// !!! Busca o key.json no diretorio do executavel (bin)
			//get { return Application.StartupPath + "\\your file name.json"; }
			get 
            {
				if (RunningPlatform() == Platform.Windows)
					FileSlash = "\\";
				else
					FileSlash = "/";

				return Application.StartupPath + FileSlash + "key.json"; 
            }

		}

        private GoogleCredential _credential;
        private GoogleCredential CreateCredential()
        {
            if (_credential != null) return _credential;
            using (var stream = new FileStream(JsonKeypath, FileMode.Open, FileAccess.Read))
            {
                string[] scopes = { VisionService.Scope.CloudPlatform };
                var credential = GoogleCredential.FromStream(stream);
                credential = credential.CreateScoped(scopes);
                _credential = credential;
                return credential;
            }
        }

        private VisionService CreateService(GoogleCredential credential)
        {
            var service = new VisionService(new BaseClientService.Initializer()
            {
                HttpClientInitializer = credential,
                ApplicationName = ApplicationName,
                GZipEnabled = true,
            });
            return service;
        }


        public void GetText(string imgPath, string language) // !!!
        {
            TextResult = JsonResult = "";
            byte[] file = File.ReadAllBytes(imgPath);

            // !!! O Path da credencial key.json deve ser o do executavel (bin)
			var credential = CreateCredential();
            var service = CreateService(credential);

			service.HttpClient.Timeout = new TimeSpan(1, 1, 1);

			BatchAnnotateImagesRequest batchRequest = new BatchAnnotateImagesRequest();
            batchRequest.Requests = new List<AnnotateImageRequest>();
            batchRequest.Requests.Add(new AnnotateImageRequest()
            {
                Features = new List<Feature>() { new Feature() { Type = "TEXT_DETECTION", MaxResults = 1 }, },
                ImageContext = new ImageContext() { LanguageHints = new List<string>() { language } },
                Image = new Image() { Content = Convert.ToBase64String(file) }
            });

			var annotate = service.Images.Annotate(batchRequest);
            BatchAnnotateImagesResponse batchAnnotateImagesResponse = annotate.Execute();
            if (batchAnnotateImagesResponse.Responses.Any())
            {
                AnnotateImageResponse annotateImageResponse = batchAnnotateImagesResponse.Responses[0];
                if (annotateImageResponse.Error != null)
                {
                    if (annotateImageResponse.Error.Message != null)
                        Error = annotateImageResponse.Error.Message;
                }
                else
                {
                    if (annotateImageResponse.TextAnnotations != null && annotateImageResponse.TextAnnotations.Any())
                    {
                        TextResult = annotateImageResponse.TextAnnotations[0].Description.Replace("\n", " "); // !!!
                        // TextResult = annotateImageResponse.TextAnnotations[0].Description.Replace("\n", "\r\n");
                        JsonResult = JsonConvert.SerializeObject(annotateImageResponse.TextAnnotations[0]);

                        // Teste de inversao de pagina
						EntityAnnotation ea = new EntityAnnotation();
                        ea = annotateImageResponse.TextAnnotations[1]; // !!!
						// ea = annotateImageResponse.TextAnnotations[0];

						List<Vertex> vertexList = ea.BoundingPoly.Vertices.ToList();
						// Calcula centro
						float centerX = 0, centerY = 0;
						for (int i = 0; i < 4; i++)
						{
                            centerX += vertexList[i].X.Value;
							centerY += vertexList[i].Y.Value;
						}
						centerX /= 4;
						centerY /= 4;

						int x0 = vertexList[0].X.Value;
						int y0 = vertexList[0].Y.Value;

                        Page_angle = 0; // Pagina normal (0 graus de rotacao)

                        if (x0 < centerX)
						{
							if (y0 < centerY)
							{
								//       0 -------- 1
								//       |          |
								//       3 -------- 2
                                Page_angle = 0; // 1
							}
							else
							{
								//       1 -------- 2
								//       |          |
								//       0 -------- 3
                                Page_angle = 270; // 6
							}
						}
						else
						{
							if (y0 < centerY)
							{
								//       3 -------- 0
								//       |          |
								//       2 -------- 1
                                Page_angle = 90; // 8
							}
							else
							{
								//       2 -------- 3
								//       |          |
								//       1 -------- 0
                                Page_angle = 180; // 3
							}
						}

                        /*
                        if (Page_angle != 0)
                            MessageBox.Show("Pagina rotacionada: " + Page_angle + " graus");
                        */
						// MessageBox.Show("(" + x1 + "," + y1 + ")" + "\r\n" + "(" + x2 + "," + y2 + ")" + "\r\n" + "(" + x3 + "," + y3 + ")" + "\r\n" + "(" + x4 + "," + y4 + ")" + "\r\n"); 
                    }
                    return;

                }

            }

            return;
            //return TextResult;
        }

    }
}
