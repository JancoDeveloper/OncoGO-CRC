# Seguridad

OncoGO-CRC es un repositorio público. No guarda historias clínicas, credenciales
ni claves de servicios.

## Cómo reportar

Si encuentra una clave, un token, una credencial o un dato personal que no
debería estar aquí, escriba a quien mantiene el repositorio por GitHub
(https://github.com/JancoDeveloper/OncoGO-CRC) o abra un aviso privado si la
plataforma se lo ofrece. No pegue el secreto en un issue público, en un commit
ni en una pull request. Si puede, indique la ruta y el commit, y el tipo de
secreto, sin el valor.

## Qué no se versiona

`.gitignore` deja fuera `.env`, variantes `.env.*`, archivos de credenciales
y de cuentas de servicio, `google-services.json`, los registros
`firebase-debug.log` y `ui-debug.log`, la carpeta `.firebase/` y los archivos
AppleDouble (`._*`) que crea macOS en discos externos.

El sitio en `web/` es estático. `firebase.json` solo dice qué carpeta publicar.
No hay claves en ese archivo ni en el sitio.

## Historial

Una clave puede seguir en commits antiguos aunque ya no esté en `main`. Este
repositorio no reescribe historia. Si eso ocurre, la medida es restringir o
rotar la clave en el servicio, no borrar commits publicados.
