IMAGE_NAME = eh-bw2:latest
PREFIX1 = weisburd
PREFIX2 = gcr.io/bw2-rare-disease
all: build push update_sha

build:
	cd .. && docker build -f docker/Dockerfile -t $(PREFIX1)/$(IMAGE_NAME) -t $(PREFIX2)/$(IMAGE_NAME) .

push:
	docker push $(PREFIX1)/$(IMAGE_NAME)
	docker push $(PREFIX2)/$(IMAGE_NAME)

update_sha:
	docker pull $(PREFIX1)/$(IMAGE_NAME) 2>&1 | grep Digest | cut -c 9- > sha256.txt
