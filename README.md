# Bitmap image manipulation

## About

This is a project I made for my Procedural Programming (C lang) class at the university.
I had to do a program that can encrypt a BMP file and decrypt it using a given secret key.

The encryption part is detailed in the `Proiect programare procedurala.pdf` file.

## The encryption part

- Read a BMP file in the intern memory
- Swap pixels order with XORSHIFT32 algorithm
- Encrypt pixels with the given key
- Save the encrypted image into BMP format.
- Decrypt the resulting image to verify that you made the encyption in a correct manner.

## The template matching part

- Read a BMP file in the intern memory
- Read the templates of the given digits (BMP files)
- Implement the template matching algorithm
- Create colored rectangles over the original image to evidentiate the patterns found.
- Eliminate the overlapping rectangles that obstruct the other rectangles already formed.

## Test yourself

1. Install CodeBlocks IDE
2. Clone this repository and open `Proiect_PP.cbp` with CodeBlocks
3. You are ready to run this project.

### Optional 

4. You can modify `cripteaza.txt` to look like: 

```txt
./input/your_image.bmp  -> Input path
./output/your_image_encrypted.bmp  -> Output path
your_secret_code.txt  -> Secret key
```

Where `./input/your_image.bmp` can be any BMP format image path and `your_secret_code.txt` is a .txt file which contains your secret encyption code.
The output image will be stored here `./output/your_image_encrypted.bmp`

5. You can modify `decripteaza.txt` to look like:

```txt
./output/output_encrypted.bmp  -> Encrypted input image
./output/output_decrypted.bmp  -> Decrypted output image
your_secret_key.txt  -> Secret key
```

Where `output_encrypted.bmp` is path to an encrypted image from the above step and `your_secret_code.txt` is a .txt file which contains your secret encyption code.
The output image will be the decrypted version of the input image and it will be stored in `./output/output_decrypted.bmp`
